/*
This tool is part of the WhiteboxTools geospatial analysis library.
Authors: Dr. John Lindsay
Created: 28/06/2017
Last Modified: 24/11/2019
License: MIT
*/

use whitebox_raster::*;
use whitebox_common::structures::Array2D;
use crate::tools::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::VecDeque;
use std::env;
use std::f64;
use std::i32;
use std::io::{Error, ErrorKind};
use std::path;

/// This tool can be used to remove the depressions in a digital elevation model (DEM), a
/// common requirement of spatial hydrological operations such as flow accumulation
/// and watershed modelling. The tool based on the efficient hybrid depression
/// breaching algorithm described by Lindsay (2016). It uses a breach-first, fill-second
/// approach to resolving continuous flowpaths through depressions.
///
/// The user may optionally (`--flat_increment`) override the default value applied to increment elevations on
/// flat areas (often formed by the subsequent depression filling operation). The default value is
/// dependent upon the elevation range in the input DEM and is generally a very small elevation value (e.g.
/// 0.001). It may be necessary to override the default elevation increment value in landscapes where there
/// are extensive flat areas resulting from depression filling (and along breach channels). Values in the range
/// 0.00001 to 0.01 are generally appropriate. Increment values that are too large can result in obvious artifacts
/// along flattened sites, which may extend beyond the flats, and values that are too small (i.e. smaller than the
/// numerical precision) may result in the presence of grid cells with no downslope neighbour in the
/// output DEM. If a flat increment value isn't specified, the output DEM will use 64-bit floating point values in order
/// to make sure that the very small elevation increment value determined will be accurately stored. Consequently,
/// it may double the storage requirements as DEMs are often stored with 32-bit precision. However, if a flat increment
/// value is specified, the output DEM will keep the same data type as the input assuming the user chose its value wisely.
///
/// In comparison with the `BreachDepressionsLeastCost` tool, this breaching method often provides a less
/// satisfactory, higher impact, breaching solution and is often less efficient. **It has been provided to users for
/// legacy reasons and it is advisable that users try the `BreachDepressionsLeastCost` tool to remove depressions from
/// their DEMs first**. The `BreachDepressionsLeastCost` tool is particularly
/// well suited to breaching through road embankments. Nonetheless, there are applications for which full depression filling
/// using the  `FillDepressions` tool may be preferred.
///
/// # Reference
/// Lindsay JB. 2016. *Efficient hybrid breaching-filling sink removal methods for
/// flow path enforcement in digital elevation models.* **Hydrological Processes**,
/// 30(6): 846–857. DOI: 10.1002/hyp.10648
///
/// # See Also
/// `BreachDepressionsLeastCost`, `FillDepressions`, `FillSingleCellPits`
pub struct BreachDepressions {
    name: String,
    description: String,
    toolbox: String,
    parameters: Vec<ToolParameter>,
    example_usage: String,
}

impl BreachDepressions {
    pub fn new() -> BreachDepressions {
        // public constructor
        let name = "BreachDepressions".to_string();
        let toolbox = "Hydrological Analysis".to_string();
        let description = "Breaches all of the depressions in a DEM using Lindsay's (2016) algorithm. This should be preferred over depression filling in most cases.".to_string();

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
            name: "Sink File".to_owned(),
            flags: vec!["--sink".to_owned()],
            description: "Sink raster file.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Raster),
            default_value: None,
            optional: false,
        });

        parameters.push(ToolParameter {
            name: "Flat increment value (z units)".to_owned(),
            flags: vec!["--flat_increment".to_owned()],
            description: "Optional elevation increment applied to flat areas.".to_owned(),
            parameter_type: ParameterType::Float,
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
            ">>.*{0} -r={1} -v --wd=\"*path*to*data*\" --dem=DEM.tif -o=output.tif",
            short_exe, name
        )
        .replace("*", &sep);

        BreachDepressions {
            name: name,
            description: description,
            toolbox: toolbox,
            parameters: parameters,
            example_usage: usage,
        }
    }
}

impl WhiteboxTool for BreachDepressions {
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
        let mut sink_file = String::new();
        let mut flat_increment = f64::NAN;

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
            } else if flag_val == "-sink" {
                sink_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
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
                }
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
            println!("Reading data...")
        };

        let mut input = Raster::new(&input_file, "r")?;


        let use_sink = !sink_file.is_empty();
        if use_sink {
            if !sink_file.contains(&sep) && !sink_file.contains("/") {
                sink_file = format!("{}{}", working_directory, sink_file);
            }
        }

        let sink: Array2D<f64> = match use_sink {
            false => Array2D::new(1, 1, -9999_f64, -9999_f64)?,
            true => {
                // if verbose { println!("Reading watershed data...") };
                let r = Raster::new(&sink_file, "r")?;
                if r.configs.rows != input.configs.rows as usize || r.configs.columns != input.configs.columns as usize {
                    return Err(Error::new(ErrorKind::InvalidInput,
                                        "The input files must have the same number of rows and columns and spatial extent."));
                }
                r.get_data_as_array2d()
            }
        };


        let start = Instant::now();
        let rows = input.configs.rows as isize;
        let columns = input.configs.columns as isize;
        let num_cells = rows * columns;
        let nodata = input.configs.nodata;
        let resx = input.configs.resolution_x;
        let resy = input.configs.resolution_y;
        let diagres = (resx * resx + resy * resy).sqrt();

        let mut output = Raster::initialize_using_file(&output_file, &input);
        let background_val = (i32::min_value() + 1) as f64;
        output.reinitialize_values(background_val);

        let small_num = if !flat_increment.is_nan() || flat_increment == 0f64 {
            output.configs.data_type = input.configs.data_type; // Assume the user knows what he's doing
            flat_increment
        } else {
            output.configs.data_type = DataType::F64; // Don't take any chances and promote to 64-bit
            let elev_digits = (input.configs.maximum as i64).to_string().len();
            let elev_multiplier = 10.0_f64.powi((6 - elev_digits) as i32);
            1.0_f64 / elev_multiplier as f64 * diagres.ceil()
        };

        let mut z: f64;
        let mut z_n: f64;
        let dx = [1, 1, 1, 0, -1, -1, -1, 0];
        let dy = [-1, 0, 1, 1, 1, 0, -1, -1];


        let mut flow_dir: Array2D<i8> = Array2D::new(rows, columns, -1, -1)?;

        /*
        Find the data edges. This is complicated by the fact that DEMs frequently
        have nodata edges, whereby the DEM does not occupy the full extent of
        the raster. One approach to doing this would be simply to scan the
        raster, looking for cells that neighbour nodata values. However, this
        assumes that there are no interior nodata holes in the dataset. Instead,
        the approach used here is to perform a region-growing operation, looking
        for nodata values along the raster's edges.
        */

        let mut queue: VecDeque<(isize, isize)> =
            VecDeque::with_capacity((rows * columns) as usize);
        for row in 0..rows {
            /*
            Note that this is only possible because Whitebox rasters
            allow you to address cells beyond the raster extent but
            return the nodata value for these regions.
            */
            queue.push_back((row, -1));
            queue.push_back((row, columns));
        }

        for col in 0..columns {
            queue.push_back((-1, col));
            queue.push_back((rows, col));
        }

        /*
        minheap is the priority queue. Note that I've tested using integer-based
        priority values, by multiplying the elevations, but this didn't result
        in a significant performance gain over the use of f64s.
        */
        let mut minheap = BinaryHeap::with_capacity((rows * columns) as usize);
        let mut num_solved_cells = 0;
        let mut zin_n: f64; // value of neighbour of row, col in input raster
        let mut zout: f64; // value of row, col in output raster
        let mut zout_n: f64; // value of neighbour of row, col in output raster
        let (mut row, mut col): (isize, isize);
        let (mut row_n, mut col_n): (isize, isize);
        while !queue.is_empty() {
            let cell = queue.pop_front().unwrap();
            row = cell.0;
            col = cell.1;
            for n in 0..8 {
                row_n = row + dy[n];
                col_n = col + dx[n];
                zin_n = input.get_value(row_n, col_n);
                zout_n = output.get_value(row_n, col_n);
                if zout_n == background_val {
                    if zin_n == nodata {
                        output.set_value(row_n, col_n, nodata);
                        queue.push_back((row_n, col_n));
                    } else {
                        output.set_value(row_n, col_n, zin_n);
                        // Push it onto the priority queue for the priority flood operation
                        minheap.push(GridCell {
                            row: row_n,
                            column: col_n,
                            priority: zin_n,
                        });
                    }
                    num_solved_cells += 1;
                }
            }

            if verbose {
                progress = (100.0_f64 * num_solved_cells as f64 / (num_cells - 1) as f64) as usize;
                if progress != old_progress {
                    println!("Progress: {}%", progress);
                    old_progress = progress;
                }
            }
        }

        // Perform the priority flood operation.
        let back_link = [4i8, 5i8, 6i8, 7i8, 0i8, 1i8, 2i8, 3i8];
        let (mut x, mut y): (isize, isize);
        let mut z_target: f64;
        let mut dir: i8;
        let mut flag: bool;


        while !minheap.is_empty() {
            let cell = minheap.pop().expect("Error during pop operation.");
            row = cell.row;
            col = cell.column;
            zout = output.get_value(row, col);
            for n in 0..8 {
                row_n = row + dy[n];
                col_n = col + dx[n];
                zout_n = output.get_value(row_n, col_n);
                if zout_n == background_val {
                    zin_n = input.get_value(row_n, col_n);
                    if zin_n != nodata {
                        flow_dir.set_value(row_n, col_n, back_link[n]);
                        output.set_value(row_n, col_n, zin_n);
                        minheap.push(GridCell {
                            row: row_n,
                            column: col_n,
                            priority: zin_n,
                        });
                        if zin_n < (zout + small_num) {
                            // Trace the flowpath back to a lower cell, if it exists.
                            x = col_n;
                            y = row_n;
                            z_target = output.get_value(row_n, col_n);
                            dir = flow_dir[(y, x)];
                            while dir >= 0 {
                                y += dy[dir as usize];
                                x += dx[dir as usize];
                                z_target -= small_num;
                                if output.get_value(y, x) > z_target {
                                    output.set_value(y, x, z_target);
                                } else {
                                    break;
                                }
                                dir = flow_dir[(y, x)];
                            }
                        }
                    } else {
                        // Interior nodata cells are still treated as nodata and are not filled.
                        output.set_value(row_n, col_n, nodata);
                        num_solved_cells += 1;
                        // region growing operation to find all attached nodata cells
                        queue.push_back((row_n, col_n));
                        while !queue.is_empty() {
                            let cell = queue.pop_front().unwrap();
                            for n2 in 0..8 {
                                let row2 = cell.0 + dy[n2];
                                let col2 = cell.1 + dx[n2];
                                if input.get_value(row2, col2) == nodata
                                    && output.get_value(row2, col2) == background_val
                                {
                                    if row2 >= 0 && row2 < rows && col2 >= 0 && col2 < columns {
                                        output.set_value(row2, col2, nodata);
                                        num_solved_cells += 1;
                                        queue.push_back((row2, col2));
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if verbose {
                num_solved_cells += 1;
                progress =
                    (100.0_f64 * num_solved_cells as f64 / (num_cells - 1) as f64) as usize;
                if progress != old_progress {
                    println!("Progress: {}%", progress);
                    old_progress = progress;
                }
            }
        }

        let elapsed_time = get_formatted_elapsed_time(start);
        output.configs.display_min = input.configs.display_min;
        output.configs.display_max = input.configs.display_max;
        output.add_metadata_entry(format!(
            "Created by whitebox_tools\' {} tool",
            self.get_tool_name()
        ));
        output.add_metadata_entry(format!("Input file: {}", input_file));
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
    // priority: usize,
    priority: f64,
}

impl Eq for GridCell {}

impl PartialOrd for GridCell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // Some(other.priority.cmp(&self.priority))
        other.priority.partial_cmp(&self.priority)
    }
}

impl Ord for GridCell {
    fn cmp(&self, other: &GridCell) -> Ordering {
        // other.priority.cmp(&self.priority)
        let ord = self.partial_cmp(other).unwrap();
        match ord {
            Ordering::Greater => Ordering::Less,
            Ordering::Less => Ordering::Greater,
            Ordering::Equal => ord,
        }
    }
}
