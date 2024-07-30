/*
This tool is part of the WhiteboxTools geospatial analysis library.
Authors: Dr. John Lindsay
Created: 01/11/2019
Last Modified: 24/11/2019
License: MIT
*/

use whitebox_raster::*;
use whitebox_common::structures::Array2D;
use crate::tools::*;
use std::cmp::Ordering;
use std::cmp::Ordering::Equal;
use std::collections::{BinaryHeap, VecDeque};
use std::env;
use std::f64;
use std::i32;
use std::io::{Error, ErrorKind};
use std::path;
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;

/// This tool can be used to perform a type of optimal depression breaching to prepare a
/// digital elevation model (DEM) for hydrological analysis. Depression breaching is a common
/// alternative to depression filling (`FillDepressions`) and often offers a lower-impact
/// solution to the removal of topographic depressions. This tool implements a method that is
/// loosely based on the algorithm described by Lindsay and Dhun (2015), furthering the earlier
/// algorithm with efficiency optimizations and other significant enhancements. The approach uses a least-cost
/// path analysis to identify the breach channel that connects pit cells (i.e. grid cells for
/// which there is no lower neighbour) to some distant lower cell. Prior to breaching and in order
/// to minimize the depth of breach channels, all pit cells are rised to the elevation of the lowest
/// neighbour minus a small heigh value. Here, the cost of a breach path is determined by the amount
/// of elevation lowering needed to cut the breach channel through the surrounding topography.
///
/// The user must specify the name of the input DEM file (`--dem`), the output breached DEM
/// file (`--output`), the maximum search window radius (`--dist`), the optional maximum breach
/// cost (`--max_cost`), and an optional flat height increment value (`--flat_increment`). Notice that **if the
/// `--flat_increment` parameter is not specified, the small number used to ensure flow across flats will be
/// calculated automatically, which should be preferred in most applications** of the tool.
/// The tool operates by performing a least-cost path analysis for each pit cell, radiating outward
/// until the operation identifies a potential breach destination cell or reaches the maximum breach length parameter.
/// If a value is specified for the optional `--max_cost` parameter, then least-cost breach paths that would require
/// digging a channel that is more costly than this value will be left unbreached. The flat increment value is used
/// to ensure that there is a monotonically descending path along breach channels to satisfy the necessary
/// condition of a downslope gradient for flowpath modelling. It is best for this value to be a small
/// value. If left unspecified, the tool with determine an appropriate value based on the range of
/// elevation values in the input DEM, **which should be the case in most applications**. Notice that the need to specify these very small elevation
/// increment values is one of the reasons why the output DEM will always be of a 64-bit floating-point
/// data type, which will often double the storage requirements of a DEM (DEMs are often store with 32-bit
/// precision). Lastly, the user may optionally choose to apply depression filling (`--fill`) on any depressions
/// that remain unresolved by the earlier depression breaching operation. This filling step uses an efficient
/// filling method based on flooding depressions from their pit cells until outlets are identified and then
/// raising the elevations of flooded cells back and away from the outlets.
///
/// The tool can be run in two modes, based on whether the `--min_dist` is specified. If the `--min_dist` flag
/// is specified, the accumulated cost (accum<sub>2</sub>) of breaching from *cell1* to *cell2* along a channel
/// issuing from *pit* is calculated using the traditional cost-distance function:
///
/// > cost<sub>1</sub> = z<sub>1</sub> - (z<sub>pit</sub> + *l* &times; *s*)
/// >
/// > cost<sub>2</sub> = z<sub>2</sub> - [z<sub>pit</sub> + (*l* + 1)*s*]
/// >
/// > accum<sub>2</sub> = accum<sub>1</sub> + *g*(cost<sub>1</sub> + cost<sub>2</sub>) / 2.0
///
/// where cost<sub>1</sub> and cost<sub>2</sub> are the costs associated with moving through *cell1* and *cell2*
/// respectively, z<sub>1</sub> and z<sub>2</sub> are the elevations of the two cells, z<sub>pit</sub> is the elevation
/// of the pit cell, *l* is the length of the breach channel to *cell1*, *g* is the grid cell distance between
/// cells (accounting for diagonal distances), and *s* is the small number used to ensure flow
/// across flats. If the `--min_dist` flag is not present, the accumulated cost is calculated as:
///
/// > accum<sub>2</sub> = accum<sub>1</sub> + cost<sub>2</sub>
///
/// That is, without the `--min_dist` flag, the tool works to minimize elevation changes to the DEM caused by
/// breaching, without considering the distance of breach channels. Notice that the value `--max_cost`, if
/// specified, should account for this difference in the way cost/cost-distances are calculated. The first cell
/// in the least-cost accumulation operation that is identified for which cost<sub>2</sub> <= 0.0 is the target
/// cell to which the breach channel will connect the pit along the least-cost path.
///
/// In comparison with the `BreachDepressions` tool, this breaching method often provides a more
/// satisfactory, lower impact, breaching solution and is often more efficient. It is therefore advisable that users
/// try the `BreachDepressionsLeastCost` tool to remove depressions from their DEMs first. This tool is particularly
/// well suited to breaching through road embankments. There are instances when a breaching solution is inappropriate, e.g.
/// when a very deep depression such as an open-pit mine occurs in the DEM and long, deep breach paths are created. Often
/// restricting breaching with the `--max_cost` parameter, combined with subsequent depression filling (`--fill`) can
/// provide an adequate solution in these cases. Nonetheless, there are applications for which full depression filling
/// using the  `FillDepressions` tool may be preferred.
///
/// # Reference
/// Lindsay J, Dhun K. 2015. Modelling surface drainage patterns in altered landscapes using LiDAR.
/// *International Journal of Geographical Information Science*, 29: 1-15. DOI: 10.1080/13658816.2014.975715
///
/// # See Also
/// `BreachDepressions`, `FillDepressions`, `CostPathway`
pub struct BreachDepressionsLeastCost {
    name: String,
    description: String,
    toolbox: String,
    parameters: Vec<ToolParameter>,
    example_usage: String,
}

impl BreachDepressionsLeastCost {
    pub fn new() -> BreachDepressionsLeastCost {
        // public constructor
        let name = "BreachDepressionsLeastCost".to_string();
        let toolbox = "Hydrological Analysis".to_string();
        let description =
            "Breaches the depressions in a DEM using a least-cost pathway method.".to_string();

        let mut parameters = vec![];
        parameters.push(ToolParameter {
            name: "Input DEM File".to_owned(),
            flags: vec!["-i".to_owned(), "--dem".to_owned()],
            description: "Input raster DEM file.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Raster),
            default_value: None,
            optional: false,
        });

        parameters.push(ToolParameter {
            name: "Output File".to_owned(),
            flags: vec!["-o".to_owned(), "--output".to_owned()],
            description: "Output raster file.".to_owned(),
            parameter_type: ParameterType::NewFile(ParameterFileType::Raster),
            default_value: None,
            optional: false,
        });

        parameters.push(ToolParameter {
            name: "Maximum Search Distance (cells)".to_owned(),
            flags: vec!["--dist".to_owned()],
            description: "Maximum search distance for breach paths in cells.".to_owned(),
            parameter_type: ParameterType::Integer,
            default_value: None,
            optional: false,
        });

        parameters.push(ToolParameter {
            name: "Minimize breach distances?".to_owned(),
            flags: vec!["--min_dist".to_owned()],
            description: "Optional flag indicating whether to minimize breach distances."
                .to_owned(),
            parameter_type: ParameterType::Boolean,
            default_value: Some("true".to_string()),
            optional: true,
        });

        parameters.push(ToolParameter {
            name: "Flat increment value (z units)".to_owned(),
            flags: vec!["--flat_increment".to_owned()],
            description: "Optional elevation increment applied to flat areas.".to_owned(),
            parameter_type: ParameterType::Float,
            default_value: None,
            optional: true,
        });

        parameters.push(ToolParameter {
            name: "Masking File".to_owned(),
            flags: vec!["--mask".to_owned()],
            description: "Masking file which indicate cells to leave untouched.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Raster),
            default_value: None,
            optional: true,
        });

        let sep: String = path::MAIN_SEPARATOR.to_string();
        let e = format!("{}", env::current_exe().unwrap().display());
        let mut parent = env::current_exe().unwrap();
        parent.pop();
        let p = format!("{}", parent.display());
        let mut short_exe = e
            .replace(&p, "")
            .replace(".exe", "")
            .replace(".", "")
            .replace(&sep, "");
        if e.contains(".exe") {
            short_exe += ".exe";
        }
        let usage = format!(
            ">>.*{0} -r={1} -v --wd=\"*path*to*data*\" --dem=DEM.tif -o=output.tif --dist=1000",
            short_exe, name
        )
        .replace("*", &sep);

        BreachDepressionsLeastCost {
            name: name,
            description: description,
            toolbox: toolbox,
            parameters: parameters,
            example_usage: usage,
        }
    }
}

impl WhiteboxTool for BreachDepressionsLeastCost {
    fn get_source_file(&self) -> String {
        String::from(file!())
    }

    fn get_tool_name(&self) -> String {
        self.name.clone()
    }

    fn get_tool_description(&self) -> String {
        self.description.clone()
    }

    fn get_tool_parameters(&self) -> String {
        match serde_json::to_string(&self.parameters) {
            Ok(json_str) => return format!("{{\"parameters\":{}}}", json_str),
            Err(err) => return format!("{:?}", err),
        }
    }

    fn get_example_usage(&self) -> String {
        self.example_usage.clone()
    }

    fn get_toolbox(&self) -> String {
        self.toolbox.clone()
    }

    fn run<'a>(
        &self,
        args: Vec<String>,
        working_directory: &'a str,
        verbose: bool,
    ) -> Result<(), Error> {
        let mut input_file = String::new();
        let mut output_file = String::new();
        let mut max_dist = 20isize;
        let mut flat_increment = f64::NAN;
        let mut mask_file = String::new();

        if args.len() == 0 {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "Tool run with no parameters.",
            ));
        }
        for i in 0..args.len() {
            let mut arg = args[i].replace("\"", "");
            arg = arg.replace("\'", "");
            let cmd = arg.split("="); // in case an equals sign was used
            let vec = cmd.collect::<Vec<&str>>();
            let mut keyval = false;
            if vec.len() > 1 {
                keyval = true;
            }
            let flag_val = vec[0].to_lowercase().replace("--", "-");
            if flag_val == "-i" || flag_val == "-input" || flag_val == "-dem" {
                input_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            } else if flag_val == "-o" || flag_val == "-output" {
                output_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            } else if flag_val == "-dist" {
                max_dist = if keyval {
                    vec[1].to_string().parse::<isize>().unwrap()
                } else {
                    args[i + 1].to_string().parse::<isize>().unwrap()
                };
            } else if flag_val == "-flat_increment" {
                flat_increment = if keyval {
                    vec[1]
                        .to_string()
                        .parse::<f64>()
                        .expect(&format!("Error parsing {}", flag_val))
                } else {
                    args[i + 1]
                        .to_string()
                        .parse::<f64>()
                        .expect(&format!("Error parsing {}", flag_val))
                };
            }  else if flag_val == "-mask" {
                mask_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            }
        }

        if verbose {
            let tool_name = self.get_tool_name();
            let welcome_len = format!("* Welcome to {} *", tool_name).len().max(28); 
            // 28 = length of the 'Powered by' by statement.
            println!("{}", "*".repeat(welcome_len));
            println!("* Welcome to {} {}*", tool_name, " ".repeat(welcome_len - 15 - tool_name.len()));
            println!("* Powered by WhiteboxTools {}*", " ".repeat(welcome_len - 28));
            println!("* www.whiteboxgeo.com {}*", " ".repeat(welcome_len - 23));
            println!("{}", "*".repeat(welcome_len));
        }

        let sep: String = path::MAIN_SEPARATOR.to_string();

        let mut progress: usize;
        let mut old_progress: usize = 1;

        if !input_file.contains(&sep) && !input_file.contains("/") {
            input_file = format!("{}{}", working_directory, input_file);
        }
        if !output_file.contains(&sep) && !output_file.contains("/") {
            output_file = format!("{}{}", working_directory, output_file);
        }

        if verbose {
            println!("Reading DEM")
        };

        let input = Arc::new(Raster::new(&input_file, "r").expect("Error reading input raster"));

        let start = Instant::now();

        let rows = input.configs.rows as isize;
        let columns = input.configs.columns as isize;
        let nodata = input.configs.nodata;
        let (mut col, mut row): (isize, isize);
        let (mut rn, mut cn): (isize, isize);
        let mut accum: f64;
        let (mut z, mut zn, mut zout): (f64, f64, f64);
        let dx = [1, 1, 1, 0, -1, -1, -1, 0];
        let dy = [-1, 0, 1, 1, 1, 0, -1, -1];
        let mut flag: bool;
        let mut num_solved: usize;
        let resx = input.configs.resolution_x;
        let resy = input.configs.resolution_y;
        let diagres = (resx * resx + resy * resy).sqrt();
        let cost_dist = [diagres, resx, diagres, resy, diagres, resx, diagres, resy];
        let mut cost1: f64;
        let mut cost2: f64;
        let mut new_cost: f64;
        let mut length: i16;
        let mut length_n: i16;
        let mut b: usize;
        let mut num_procs = num_cpus::get() as isize;
        let configs = whitebox_common::configs::get_configs()?;
        let max_procs = configs.max_procs;
        if max_procs > 0 && max_procs < num_procs {
            num_procs = max_procs;
        }


        let mut mask = Array2D::new(rows, columns, 0i8, 0i8)?;
        if !mask_file.is_empty() {
            println!("Reading mask");
            let mask_f64 = Raster::new(&mask_file, "r").expect("Error reading mask raster file");
        
            // Make sure the mask file has the same size of the input file (resolution is assumed to be the same)
        if input.configs.rows != mask_f64.configs.rows  || input.configs.columns != mask_f64.configs.columns {
            return Err(Error::new(ErrorKind::InvalidInput,
                                "The mask and dem raster files must have the same number of rows and columns and spatial extent."));
        }

        // Transposition des valeurs dans un masque 8bit
        // Tout ce qui n'est pas 0 fera partie du masque
        for row in 0..rows {
            for col in 0..columns {
                z = mask_f64.get_value(row, col);
                if z != 0.0 {
                        mask.set_value(row, col, 1i8);
                }
            }
        }
        
        drop(mask_f64);
        }
        let mask = Arc::new(mask.clone());


        // 1 - Procéder à l'ajustement des pixels en bordure du masque pour éviter,
        //     On parle ici d'un abaissement ou élévation des pixels pour éviter tout débordement ou attirance
        //     des pixels sous le masque. Penser que je dois respecter le sens d'écoulement d'origine des pixels
        //     à l'intérieur du masque (cas de l'exutoire par exemple). je dois donc considérer le cas où l'eau
        //     entre naturellement dans la zone d'intérêt







        // Grossi d'un pixel le masque afin d'englober les pixels de bordure du masque, soit les pixels qui ne
        // doivent pas être modifiés par le remplissage ou le bréchage afin de ne pas compromettre l'hydrocohérence
        // avec les bassins versants voisins
        let mut mask_expand = Array2D::new(rows, columns, 0i8, 0i8)?;
        if !mask_file.is_empty() {
            println!("Expanding mask");
            let (tx, rx) = mpsc::channel();
            for tid in 0..num_procs {
                let mask = mask.clone();
                let tx = tx.clone();
                thread::spawn(move || {
                    let mut z_n: i8;
                    let dx = [1, 1, 1, 0, -1, -1, -1, 0, 0];
                    let dy = [-1, 0, 1, 1, 1, 0, -1, -1, 0];
                    for row in (0..rows).filter(|r| r % num_procs == tid) {
                        let mut data = vec![0i8; columns as usize];
                        for col in 0..columns {
                            for i in 0..9 {
                                z_n = mask.get_value(row + dy[i], col + dx[i]);
                                if z_n == 1i8 {
                                    data[col as usize] = 1i8;
                                    break;
                                }
                            }
                        }
                        tx.send((row, data)).unwrap();
                    }
                });
            }

            for _ in 0..rows {
                let (row, data) = rx.recv().expect("Error receiving data from thread.");
                mask_expand.set_row_data(row, data);
            }
        }
        let mask_expand = Arc::new(mask_expand.clone());


        let small_num = flat_increment;
        let mut output = Raster::initialize_using_file(&output_file, &input);
        // Even if the input is f32, the output will need to be 64-bit to represent the small elevation differences
        output.configs.data_type = DataType::F64;
        let display_min = input.configs.display_min;
        let display_max = input.configs.display_max;

        // Raise pit cells to minimize the depth of breach channels.
        // Only cells outside the mask will be checked and only the
        // cells that aren't on the edge of the mask will be raised.
        println!("Finding pits");
        let (tx, rx) = mpsc::channel();
        for tid in 0..num_procs {
            let input = input.clone();
            let mask = mask.clone();
            let mask_expand = mask_expand.clone();
            let tx = tx.clone();
            thread::spawn(move || {
                let (mut z, mut zn, mut min_zn, mut zm, mut zme): (f64, f64, f64, i8, i8);
                let mut is_pit: bool;
                let dx = [1, 1, 1, 0, -1, -1, -1, 0];
                let dy = [-1, 0, 1, 1, 1, 0, -1, -1];
                for row in (0..rows).filter(|r| r % num_procs == tid) {
                    let mut data = input.get_row_data(row);
                    let mut pits = vec![];
                    for col in 0..columns {
                        z = input.get_value(row, col);
                        if z != nodata {
                            // Proceed only if the cell isn't inside the mask
                            zm = mask.get_value(row, col);
                            if zm == 0i8 {
                                is_pit = true;
                            min_zn = f64::INFINITY;
                            for n in 0..8 {
                                zn = input.get_value(row + dy[n], col + dx[n]);
                                    if zn < z {
                                        // There's a lower neighbour
                                        is_pit = false;
                                        break;
                                }
                                if zn == nodata {
                                        // It's an edge cell
                                        is_pit = false;
                                    break;
                                }
                                    if zn < min_zn {
                                        min_zn = zn;
                                }
                            }
                                if is_pit {
                                    // Raise pit only if it isn't inside the expanded mask
                                    zme = mask_expand.get_value(row, col);
                                    if zme == 0i8 {
                                data[col as usize] = min_zn - small_num;
                                    }
                                pits.push((row, col, z));
                            }
                        }
                    }
                    }
                    tx.send((row, data, pits)).unwrap();
                }
            });
        }

        let mut undefined_flow_cells: Vec<(isize, isize, f64)> = vec![];
        for _ in 0..rows {
            let (row, data, mut pits) = rx.recv().expect("Error receiving data from thread.");
            output.set_row_data(row, data);
            undefined_flow_cells.append(&mut pits);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        // We need to visit and (potentially) solve each undefined-flow cell in order from lowest //
        // to highest. This is because some higher pits can be solved, or partially solved using  //
        // the breach paths of lower pits.                                                        //
        ////////////////////////////////////////////////////////////////////////////////////////////

        /* Vec is a stack and so if we want to pop the values from lowest to highest, we need to sort
        them from highest to lowest. */
        undefined_flow_cells.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Equal));
        undefined_flow_cells.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Equal));
        undefined_flow_cells.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(Equal));
        let num_deps = undefined_flow_cells.len();
        println!("Number of pits found: {}", num_deps);
        num_solved = 0;
        let backlink_dir = [4i8, 5, 6, 7, 0, 1, 2, 3];
        let mut backlink: Array2D<i8> = Array2D::new(rows, columns, -1, -2)?;
        let mut path_length: Array2D<i16> = Array2D::new(rows, columns, 0, -1)?;
        let mut scanned_cells = vec![];
        let max_length = max_dist as i16;
        let filter_size = ((max_dist * 2 + 1) * (max_dist * 2 + 1)) as usize;
        let mut minheap = BinaryHeap::with_capacity(filter_size);


        // Remplissage de la couche des pixels rencontrés avec la couche de masque (gonflé) afin d'être
        // certain que l'algorithme ne touche pas aux pixels du masque
        let mut encountered: Array2D<i8> = Array2D::new(rows, columns, 0, -1)?;
        for row in 0..rows {
            for col in 0..columns {
                encountered.set_value(row, col, mask_expand.get_value(row, col))
            }
        }



        ////////////////////////////////
        // Start of least cost breaching
        ////////////////////////////////
        while let Some(cell) = undefined_flow_cells.pop() {
            row = cell.0;
            col = cell.1;
            z = output.get_value(row, col);

            // Is it still a pit cell? It may have been solved during a previous depression solution.
            flag = true;
            for n in 0..8 {
                zn = output.get_value(row + dy[n], col + dx[n]);
                if zn < z && zn != nodata {
                    // It has a lower non-nodata cell
                    // Resolving some other pit cell resulted in a solution for this one.
                    num_solved += 1;
                    flag = false;
                    break;
                }
            }
            if flag {
                // Perform the cost-accumulation operation.
                encountered.set_value(row, col, 1i8);
                if !minheap.is_empty() {
                    minheap.clear();
                }
                minheap.push(GridCell {
                    row: row,
                    column: col,
                    priority: 0f64,
                });
                scanned_cells.push((row, col));
                flag = true;
                while !minheap.is_empty() && flag {
                    let cell2 = minheap.pop().expect("Error during pop operation.");
                    accum = cell2.priority;
                    length = path_length.get_value(cell2.row, cell2.column);
                    zn = output.get_value(cell2.row, cell2.column);
                    cost1 = zn - z + length as f64 * small_num;
                    for n in 0..8 {
                        cn = cell2.column + dx[n];
                        rn = cell2.row + dy[n];
                        if encountered.get_value(rn, cn) != 1i8 {
                            scanned_cells.push((rn, cn));
                            // not yet encountered
                            length_n = length + 1;
                            path_length.set_value(rn, cn, length_n);
                            backlink.set_value(rn, cn, backlink_dir[n]);
                            zn = output.get_value(rn, cn);
                            zout = z - (length_n as f64 * small_num);
                            if zn > zout && zn != nodata {
                                cost2 = zn - zout;
                                new_cost = accum + (cost1 + cost2) / 2f64 * cost_dist[n];  // on considère qu'on veut toujours minimiser la distance
                                encountered.set_value(rn, cn, 1i8);
                                if length_n <= max_length {
                                    minheap.push(GridCell {
                                        row: rn,
                                        column: cn,
                                        priority: new_cost,
                                    });
                                }
                            } else if zn <= zout || zn == nodata {
                                // We're at a cell that we can breach to
                                while flag {
                                    // Find which cell to go to from here
                                    if backlink.get_value(rn, cn) > -1i8 {
                                        b = backlink.get_value(rn, cn) as usize;
                                        rn += dy[b];
                                        cn += dx[b];
                                        zn = output.get_value(rn, cn);
                                        length = path_length.get_value(rn, cn);
                                        zout = z - (length as f64 * small_num);
                                        if zn > zout {
                                            output.set_value(rn, cn, zout);
                                        }
                                    } else {
                                        flag = false;
                                    }
                                }
                                num_solved += 1;
                                flag = false;
                                break; // don't check any more neighbours.
                            }
                        }
                    }
                }

                // clear the intermediate rasters
                while let Some(cell2) = scanned_cells.pop() {
                    backlink.set_value(cell2.0, cell2.1, -1i8);
                    encountered.set_value(cell2.0, cell2.1, 0i8);
                    path_length.set_value(cell2.0, cell2.1, 0i16);
                }

                // Reset to initial value of the on encountered array for the first cell
                // This is needed for cells on the edge of the mask in order to prevent the algorithm
                // editing the cells value on a subsequent breaching operation
                encountered.set_value(row, col, mask_expand.get_value(row, col));
            }

            if verbose {
                progress = (100.0_f64 * (1f64 - (undefined_flow_cells.len()) as f64 / (num_deps - 1) as f64)) as usize;
                if progress != old_progress {
                    println!("Breaching: {}%", progress);
                    old_progress = progress;
                }
            }
        }
        if verbose {
            println!("Number of pits solved: {}", num_solved);
        }


                }
            }

            if verbose {
                progress = (100.0_f64
                    * (1f64 - (undefined_flow_cells.len()) as f64 / (num_deps - 1) as f64))
                    as usize;
                if progress != old_progress {
                    println!("Breaching: {}%", progress);
                    old_progress = progress;
                }
            }
        }
        if verbose {
            println!("Num. solved pits: {}", num_solved);
            }


        let elapsed_time = get_formatted_elapsed_time(start);
        output.configs.display_min = display_min;
        output.configs.display_max = display_max;
        output.add_metadata_entry(format!(
            "Created by whitebox_tools\' {} tool",
            self.get_tool_name()
        ));
        output.add_metadata_entry(format!("Input file: {}", input_file));
        output.add_metadata_entry(format!("Maximum search distance: {}", max_dist));
        output.add_metadata_entry(format!("Flat elevation increment: {}", small_num));
        output.add_metadata_entry(format!("Elapsed Time (excluding I/O): {}", elapsed_time));

        if verbose {
            println!("Saving data...")
        };
        let _ = match output.write() {
            Ok(_) => {
                if verbose {
                    println!("Output file written")
                }
            }
            Err(e) => return Err(e),
        };
        if verbose {
            println!(
                "{}",
                &format!("Elapsed Time (excluding I/O): {}", elapsed_time)
            );
        }

        Ok(())
    }
}

#[derive(PartialEq, Debug)]
struct GridCell {
    row: isize,
    column: isize,
    priority: f64,
}

impl Eq for GridCell {}

impl PartialOrd for GridCell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.priority.partial_cmp(&self.priority)
    }
}

impl Ord for GridCell {
    fn cmp(&self, other: &GridCell) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[derive(PartialEq, Debug)]
struct GridCell2 {
    row: isize,
    column: isize,
    z: f64,
    priority: f64,
}

impl Eq for GridCell2 {}

impl PartialOrd for GridCell2 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.priority.partial_cmp(&self.priority)
    }
}

impl Ord for GridCell2 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
