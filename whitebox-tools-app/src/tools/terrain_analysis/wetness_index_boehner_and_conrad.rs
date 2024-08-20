/*
This tool is part of the WhiteboxTools geospatial analysis library.
Authors: Jean-François Bourdon
Created: 05/08/2024
Last Modified: 05/08/2024
License: MIT
*/

use whitebox_raster::*;
use whitebox_common::structures::Array2D;
use crate::tools::*;
use std::cmp::Ordering::Equal;
use num_cpus;
use std::env;
use std::f32;
use std::f64;
use std::io::{Error, ErrorKind};
use std::path;
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;


/// This tool can be used to calculate the topographic wetness index commonly used in the TOPMODEL rainfall-runoff framework.
/// The index describes the propensity for a site to be saturated to the surface given its contributing area and slope
/// characteristics. It is calculated as:
///
/// > WI = Ln(Area / tan(Slope))
///
/// Where `Area` is the catchment area (i.e. the upslope contributing area) estimated using a multiple-flow-direction
/// accumulation algorithm. The initial catchment area values are however modified iteratively according to the highest
/// contributing area of the neighbouring cells weighted by their slope. This iterative process better predicts potential
/// saturation for cells situated in valley floors with a small vertical distance to a channel compared to the standard
/// TWI calculation from `WetnessIndex` at the cost of longer computing times.
///
/// The DEM must have been hydrologically corrected to remove all spurious depressions and flat areas. DEM pre-processing
/// is usually achieved using either the `BreachDepressions` (also `BreachDepressionsLeastCost`) or `FillDepressions` tool.
/// The output raster is of the float data type and continuous data scale.
/// 
/// Derived from the C++ implementation of the *SAGA Wetness Index* tool by Olaf Conrad in SAGA GIS.
///
/// # References
/// Boehner, J., Koethe, R. Conrad, O., Gross, J., Ringeler, A., & Selige, T. 2002.
/// *Soil Regionalisation by Means of Terrain Analysis and Process Parameterisation.*
/// In: Micheli, E., Nachtergaele, F., Montanarella, L. [Ed.]: Soil Classification 2001.
/// European Soil Bureau, Research Report No. 7, EUR 20398 EN, Luxembourg: 213-222.
/// 
/// Boehner, J., & Selige, T. 2006. *Spatial prediction of soil attributes using terrain
/// analysis and climate regionalisation.* In: Boehner, J., McCloy, K.R., Strobl, J.
/// [Eds.]: SAGA - Analysis and Modelling Applications, Goettinger Geographische Abhandlungen,
/// Goettingen: 13-28.
/// 
/// See Also
/// `WetnessIndex`, `BreachDepressionsLeastCost`, `FillDepressions`
pub struct WetnessIndexBoehnerAndConrad {
    name: String,
    description: String,
    toolbox: String,
    parameters: Vec<ToolParameter>,
    example_usage: String,
}

impl WetnessIndexBoehnerAndConrad {
    pub fn new() -> WetnessIndexBoehnerAndConrad {
        // public constructor
        let name = "WetnessIndexBoehnerAndConrad".to_string();
        let toolbox = "Geomorphometric Analysis".to_string();
        let description =
            "Calculates the topographic wetness index from SAGA GIS.".to_string();

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
            name: "Weights File".to_owned(),
            flags: vec!["--weights".to_owned()], 
            description: "Weights raster file for initial catchment area.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Raster),
            default_value: None,
            optional: true,
        });

        parameters.push(ToolParameter{
            name: "Area Type".to_owned(), 
            flags: vec!["--area_type".to_owned()], 
            description: "Area type; one of 'total catchment area', 'square root of catchment area', or 'specific catchment area (default)'.".to_owned(),
            parameter_type: ParameterType::OptionList(vec!["total catchment area".to_owned(), "square root of catchment area".to_owned(), "specific catchment area".to_owned()]),
            default_value: Some("specific catchment area".to_owned()),
            optional: true
        });

        parameters.push(ToolParameter{
            name: "Slope Type".to_owned(), 
            flags: vec!["--slope_type".to_owned()], 
            description: "Slope type; one of 'local slope' or 'catchment slope (default)'.".to_owned(),
            parameter_type: ParameterType::OptionList(vec!["local slope".to_owned(), "catchment slope".to_owned()]),
            default_value: Some("catchment slope".to_owned()),
            optional: true
        });

        parameters.push(ToolParameter {
            name: "Suction".to_owned(),
            flags: vec!["--suction".to_owned()],
            description: "Optional suction factor (default is 10.0).".to_owned(),
            parameter_type: ParameterType::Float,
            default_value: Some("10".to_owned()),
            optional: true,
        });

        parameters.push(ToolParameter {
            name: "Slope weight".to_owned(),
            flags: vec!["--slope_weight".to_owned()],
            description: "Optional slope weight for index calculation (default is 1.0).".to_owned(),
            parameter_type: ParameterType::Float,
            default_value: Some("1".to_owned()),
            optional: true,
        });

        parameters.push(ToolParameter {
            name: "Minimum slope (in degrees)".to_owned(),
            flags: vec!["--slope_min".to_owned()],
            description: "Optional minimum slope for index calculation (default is 0.0).".to_owned(),
            parameter_type: ParameterType::Float,
            default_value: Some("0".to_owned()),
            optional: true,
        });

        parameters.push(ToolParameter {
            name: "Slope offset (in degrees)".to_owned(),
            flags: vec!["--slope_offset".to_owned()],
            description: "Optional slope offset for index calculation (default is 0.1).".to_owned(),
            parameter_type: ParameterType::Float,
            default_value: Some("0.1".to_owned()),
            optional: true,
        });

        parameters.push(ToolParameter {
            name: "MFD convergence".to_owned(),
            flags: vec!["--mfd_convergence".to_owned()],
            description: "Optional MFD convergence parameter (default is 1.1).".to_owned(),
            parameter_type: ParameterType::Float,
            default_value: Some("1.1".to_owned()),
            optional: true,
        });

        parameters.push(ToolParameter {
            name: "Output File".to_owned(),
            flags: vec!["-o".to_owned(), "--output".to_owned()],
            description: "Output raster file.".to_owned(),
            parameter_type: ParameterType::NewFile(ParameterFileType::Raster),
            default_value: None,
            optional: false,
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
        let usage = format!(">>.*{0} -r={1} -v --wd=\"*path*to*data*\" --dem='dem.tif' -o=output.tif", short_exe, name).replace("*", &sep);

        WetnessIndexBoehnerAndConrad {
            name: name,
            description: description,
            toolbox: toolbox,
            parameters: parameters,
            example_usage: usage,
        }
    }
}

impl WhiteboxTool for WetnessIndexBoehnerAndConrad {
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
        let mut s = String::from("{\"parameters\": [");
        for i in 0..self.parameters.len() {
            if i < self.parameters.len() - 1 {
                s.push_str(&(self.parameters[i].to_string()));
                s.push_str(",");
            } else {
                s.push_str(&(self.parameters[i].to_string()));
            }
        }
        s.push_str("]}");
        s
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
        let mut dem_file = String::new();
        let mut weights_file = String::new();
        let mut output_file = String::new();
        let mut suction_weight = 10_f32;
        let mut slope_weight = 1_f32;
        let mut area_type = 2_isize;
        let mut slope_type = 1_isize;
        let mut slope_min = 0_f32;
        let mut slope_offset = 0.1_f32;
        let mut mfd_convergence = 1.1_f64;


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
                dem_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            
            } else if flag_val == "-weights" {
                weights_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                }

            } else if flag_val == "-o" || flag_val == "-output" {
                output_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                }

            } else if flag_val == "-area_type" {
                let area_type_flag = if keyval {
                    vec[1].to_lowercase()
                } else {
                    args[i + 1].to_lowercase()
                };
                area_type = match area_type_flag.as_str() {
                    "total" => 0_isize, // total catchment area
                    "square" => 1_isize, // square root of catchment area
                    "specific" => 2_isize, // specific catchment area
                    _ => panic!("Invalid 'area_type' parameter"),
                };

            } else if flag_val == "-slope_type" {
                let slope_type_flag = if keyval {
                    vec[1].to_lowercase()
                } else {
                    args[i + 1].to_lowercase()
                };
                slope_type = match slope_type_flag.as_str() {
                    "local" => 0_isize, // local slope
                    "catchment" => 1_isize, // catchment slope
                    _ => panic!("Invalid 'slope_type' parameter"),
                };

            } else if flag_val == "-suction" {
                suction_weight = if keyval {
                    vec[1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                } else {
                    args[i + 1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                }

            } else if flag_val == "-slope_weight" {
                slope_weight = if keyval {
                    vec[1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                } else {
                    args[i + 1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                }

            } else if flag_val == "-slope_min" {
                slope_min = if keyval {
                    vec[1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                } else {
                    args[i + 1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                }

            } else if flag_val == "-slope_offset" {
                slope_offset = if keyval {
                    vec[1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                } else {
                    args[i + 1]
                        .to_string()
                        .parse::<f32>()
                        .expect(&format!("Error parsing {}", flag_val))
                }

            } else if flag_val == "-mfd_convergence" {
                mfd_convergence = if keyval {
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

        if !dem_file.contains(&sep) && !dem_file.contains("/") {
            dem_file = format!("{}{}", working_directory, dem_file);
        }

        if !output_file.contains(&sep) && !output_file.contains("/") {
            output_file = format!("{}{}", working_directory, output_file);
        }

        if verbose {
            println!("Reading DEM data...")
        };
        let dem = Arc::new(Raster::new(&dem_file, "r")?);

        let start = Instant::now();
        let rows = dem.configs.rows as isize;
        let columns = dem.configs.columns as isize;
        let nodata = dem.configs.nodata;
        let resx = dem.configs.resolution_x;
        let resy = dem.configs.resolution_y;

        let mut num_procs = num_cpus::get() as isize;
        let configs = whitebox_common::configs::get_configs()?;
        let max_procs = configs.max_procs;
        if max_procs > 0 && max_procs < num_procs {
            num_procs = max_procs;
        }
        

        let mut m_weights = Array2D::new(rows, columns, 1f32, -1f32)?;
        if !weights_file.is_empty() { 
            if verbose {
                println!("Reading weights data...")
            };
            if !weights_file.contains(&sep) && !weights_file.contains("/") {
                weights_file = format!("{}{}", working_directory, weights_file);
            }
            let r = Raster::new(&weights_file, "r")?;
            if r.configs.rows != rows as usize || r.configs.columns != columns as usize {
                return Err(Error::new(ErrorKind::InvalidInput,
                                    "The input files must have the same number of rows and columns and spatial extent."));
            }
            m_weights = r.get_data_as_f32_array2d();
        }



        ///////////
        // Calcul des accumulations de flux selon MFD
        // DEBUT DE get_area()
        ///////////
        if verbose {
            println!("Calculate initial slope and suction matrices...")
        };

        // Calcul de la matrice initiale de pente, de la matrice de suction
        // ainsi que de l'index d'élévation
        let (tx, rx) = mpsc::channel();
        for tid in 0..num_procs {
            let dem = dem.clone();
            let tx = tx.clone();
            thread::spawn(move || {
                for row in (0..rows).filter(|r| r % num_procs == tid) {
                    let mut vec_cells = Vec::<(isize, isize, f64)>::with_capacity(columns as usize);
                    let mut vec_slope = vec![-1f32; columns as usize];
                    let mut vec_suction = vec![-1f32; columns as usize];
                    for col in 0..columns {
                        let z = dem.get_value(row, col);
                        if z != nodata {
                            // Calcul de la pente initiale
                            let slope = get_gradient(&dem, row, col, z, resx);
                            vec_slope[col as usize] = slope;

                            // Calcul de la suction
                            let t_param = suction_weight.powf(slope_weight * slope);
                            vec_suction[col as usize] = (1.0 / t_param).powf(t_param.exp());

                            // Ajout à l'index d'élévation
                            vec_cells.push((row, col, z));

                        }
                    }
                    tx.send((row, vec_slope, vec_suction, vec_cells)).unwrap();
                }
            });
        }

        let mut m_suction: Array2D<f32> = Array2D::new(rows, columns, 0f32, -1f32)?;
        let mut m_slope: Array2D<f32> = Array2D::new(rows, columns, 0f32, -1f32)?;
        let mut cells_ordered = Vec::<(isize, isize, f64)>::with_capacity((rows * columns) as usize);
        for _ in 0..rows {
            let (row, vec_slope, vec_suction, mut vec_cells) = rx.recv().expect("Error receiving data from thread.");
            m_slope.set_row_data(row, vec_slope);
            m_suction.set_row_data(row, vec_suction);
            cells_ordered.append(&mut vec_cells);
        }

        // In order to pop the values from highest to lowest, we need to sort them from lowest to highest.
        // To ensure constant ordering from run to run (due to multiprocessing), values are first sorted by row and column
        if verbose {
            println!("Sorting cells...")
        };
        cells_ordered.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Equal));
        cells_ordered.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Equal));
        cells_ordered.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(Equal));



        // Calcul du MFD initial
        // N'est finalement pas parallélisable à cause de la modification
        // itérative de "m_area" et "m_slope"
        if verbose {
            println!("Calculate initial MFD...")
        };
        let mut m_area: Array2D<f32> = Array2D::new(rows, columns, -1f32, -1f32)?;
        let dcol = [0, 1, 1, 1, 0, -1, -1, -1];
        let drow = [1, 1, 0, -1, -1, -1, 0, 1];
        let diagres = (resx * resx + resy * resy).sqrt();
        let grid_lengths = [resy, diagres, resx, diagres, resy, diagres, resx, diagres];

        let nb_cells = cells_ordered.len();
        let mut ii = 0_usize;


        while let Some(cell) = cells_ordered.pop() {
            let row: isize = cell.0;
            let col: isize = cell.1;
            let z: f64 = cell.2;

            // Ajustement initial de l'accumulation (non parallélisable à cause de la modification itérative de "m_area" plus loin)
            // L'ajout de 1 est simplement pour contrebalancer les valeurs de départ de m_area qui sont de -1. Considérant que je passe
            // à travers toutes les cellules valides avec ce while, ça me permet de conserver à -1 uniquement les cellules NoData
            let area = m_area.get_value(row, col) + m_weights.get_value(row, col) + 1f32;
            m_area.set_value(row, col, area);

            // Ajustement initial de la pente du catchment en fonction de l'accumulation
            let slope = m_slope.get_value(row, col);
            m_slope.set_value(row, col, slope / area);


            // Ajustement final de l'accumulation et de la pente du catchment
            let mut dz = vec![0_f32; 8];
            let mut dz_sum = 0_f32;

            for ii in 0..8 {
                let row_n = row + drow[ii];
                let col_n = col + dcol[ii];
                let z_n = dem.get_value(row_n, col_n);
                if z_n != nodata {
                    let d = z - z_n;
                    if d > 0.0 {
                        dz[ii] = (d / grid_lengths[ii]).atan().powf(mfd_convergence) as f32;
                        dz_sum += dz[ii];
                    }
                }
            }

            if dz_sum > 0.0 {
                for ii in 0..8 {
                    if dz[ii] > 0.0 {
                        let row_n = row + drow[ii];
                        let col_n = col + dcol[ii];
                        let z_n = dem.get_value(row_n, col_n);
                        if z_n != nodata {
                            m_area.increment(row_n, col_n, area * dz[ii] / dz_sum);
                            m_slope.increment(row_n, col_n, slope * dz[ii] / dz_sum);
                        }
                    }
                }
            }


            if verbose {
                ii += 1;
                progress = (100.0_f64 * ii as f64 / (nb_cells - 1) as f64) as usize;
                if progress != old_progress {
                    println!("Initial MFD: {}%", progress);
                    old_progress = progress;
                }
            }


        }


        // Ajustement de l'accumulation en fonction de la taille de cellule
        // Pas pertinent à paralléliser, il risque d'y avoir utimement plus de visites
        // de cellules
        if verbose {
            println!("Adjust MFD to cell area...")
        };
        let cell_area = (resx * resy) as f32;
        for row in 0..rows {
            for col in 0..columns {
                let z = dem.get_value(row, col);
                if z != nodata {
                    m_area.set_value(row, col, m_area.get_value(row, col) * cell_area);
                }
            }
        }

        // FIN DE get_area()




        // Ajustement des valeurs d'accumulation en
        // fonction de la couche d'accumulation utilisateur


        // Calcul itératif de l'accumulation de flux modifiée
        if verbose {
            println!("Modify MFD...")
        };
        let m_amod = get_modified(m_area, m_suction);


        // Calcul classique du TWI
        if verbose {
            println!("Calculate topographic wetness index...")
        };
        let mut twi = Raster::initialize_using_file(&output_file, &dem);
        twi.configs.data_type = DataType::F32;
        get_twi(&mut twi, m_amod, m_slope, dem, area_type, slope_type, slope_min, slope_offset);







        
        let elapsed_time = get_formatted_elapsed_time(start);
        
        twi.add_metadata_entry(format!(
            "Created by whitebox_tools\' {} tool",
            self.get_tool_name()
        ));
        twi.add_metadata_entry(format!("Elapsed Time (excluding I/O): {}", elapsed_time));

        if verbose {
            println!("Saving data...")
        };
        let _ = match twi.write() {
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




fn get_gradient<'a>(dem: &'a Raster, row: isize, col: isize, z:f64, cellsize: f64) -> f32 {
    let nodata = dem.configs.nodata;
    let mut dz = vec![0.0_f64; 4];

    // Calcul de la pente initiale
    let dcol = [0, 1, 0,-1];
    let drow = [1, 0,-1, 0];

    for ii in 0..4 {
        let col_to = col + dcol[ii];
        let row_to = row + drow[ii];
        let col_from = col - dcol[ii];
        let row_from = row - drow[ii];

        let zn_to = dem.get_value(row_to, col_to);
        let zn_from = dem.get_value(row_from, col_from);
        if zn_to != nodata {
            dz[ii] = zn_to - z;
        } else if zn_from != nodata {
            dz[ii] = z - zn_from;
        }
    }

    let g = (dz[0] - dz[2]) / (2.0 * cellsize);
    let h = (dz[1] - dz[3]) / (2.0 * cellsize);

    let slope = (g*g + h*h).sqrt().atan() as f32;
    // aspect = (-h/-g).atan();

    return slope;
}


fn get_modified(m_area_ini: Array2D<f32>, m_suction: Array2D<f32>) -> Array2D<f32> {
    let rows = m_area_ini.rows as isize;
    let columns = m_area_ini.columns as isize;
    let num_procs = num_cpus::get() as isize;

    let mut m_amod = m_area_ini.duplicate();
    let mut m_area = m_area_ini.duplicate();

    let mut nb_changes = 1usize;
    let mut iteration = 0usize;

    while nb_changes > 0 {
        iteration += 1;
        nb_changes = 0;

        // Boucle parallélisable même si "area_mod" nécessite les valeurs voisine dans get_local_maximum()
        // car il y a ultimement convergence. Ça prend juste quelques itération supplémentaires.
        // Autre chose à tester: intégrer directement ici la fonction get_local_maximum
        // pour voir s'il y a un impact sur la performance (éviterait une déclaration répétitive
        // de variables, mais le compilateur voit peut-être les choses autrement).
        for row in 0..rows {
            for col in 0..columns {
                let area_mod = m_suction.get_value(row, col) * get_local_maximum(&m_area, row, col);
                let area = m_area.get_value(row, col);
                if area_mod > area {
                    nb_changes += 1;
                    m_area.set_value(row, col, area_mod);
                }
            }
        }

        if nb_changes > 0 {
            nb_changes = 0;

            // Boucle parallélisable
            for row in 0..rows {
                for col in 0..columns {
                    let area = m_area.get_value(row, col);
                    if area != m_amod.get_value(row, col) {
                        nb_changes += 1;
                        m_amod.set_value(row, col, area);
                    }
                }
            }
        }


        println!("pass {iteration} ({nb_changes} > 0)");
    }


    println!("post-processing...");

    let (tx, rx) = mpsc::channel();
    for tid in 0..num_procs {
        let m_area_ini = m_area_ini.clone();
        let m_area = m_area.clone();
        let dcol = [0, 1, 1, 1, 0, -1, -1, -1];
        let drow = [1, 1, 0, -1, -1, -1, 0, 1];
        let tx = tx.clone();
        thread::spawn(move || {
            for row in (0..rows).filter(|r| r % num_procs == tid) {
                let mut vec_amod = m_area.get_row_data(row);
                for col in 0..columns {
                    if m_area_ini.get_value(row, col) != m_area.nodata {
                        let mut area_modified = false;
                        let mut n = 1_isize;
                        let mut z = vec_amod[col as usize];
                        for ii in 0..8 {
                            let row_n = row + drow[ii];
                            let col_n = col + dcol[ii];
                            let area_ini = m_area_ini.get_value(row_n, col_n);
                            if area_ini != m_area.nodata {
                                let area = m_area.get_value(row_n, col_n);
                                if area > area_ini {
                                    area_modified = true;
                                }
                                n += 1;
                                z += area;
                            }
                        }
                        if area_modified {
                            vec_amod[col as usize] = z / n as f32;
                        }
                    }
                }
                tx.send((row, vec_amod)).unwrap();
            }
        });
    }

    for _ in 0..rows {
        let (row, vec_amod) = rx.recv().expect("Error receiving data from thread.");
        m_amod.set_row_data(row, vec_amod);
    }

    return m_amod;
}


fn get_local_maximum<'a>(m_grid: &'a Array2D<f32>, row: isize, col: isize) -> f32 {
    let nodata = m_grid.nodata;
    let dcol = [0, 1, 1, 1, 0, -1, -1, -1];
    let drow = [1, 1, 0, -1, -1, -1, 0, 1];
    let (mut row_n, mut col_n): (isize, isize);
    let mut z_n: f32;

    let mut val_max = m_grid.get_value(row, col);

    for ii in 0..8 {
        row_n = row + drow[ii];
        col_n = col + dcol[ii];
            z_n = m_grid.get_value(row_n, col_n);
            if z_n != nodata {
                if z_n > val_max {
                    val_max = z_n;
                }
            }
    }
    return val_max;
}


fn get_twi<'a>(twi: &'a mut Raster, m_amod: Array2D<f32>, m_slope: Array2D<f32>, dem: Arc<Raster>, area_type: isize, slope_type: isize, slope_min: f32, slope_offset: f32) {
    let rows = dem.configs.rows as isize;
    let columns = dem.configs.columns as isize;
    let nodata = dem.configs.nodata;
    let cellsize = dem.configs.resolution_x * dem.configs.resolution_y;
    let num_procs = num_cpus::get() as isize;

    let slope_min_rad = slope_min.to_radians();
    let slope_offset_rad = slope_offset.to_radians();

    let (tx, rx) = mpsc::channel();
    for tid in 0..num_procs {
        let dem = dem.clone();
        let m_slope = m_slope.clone();
        let m_amod = m_amod.clone();
        let tx = tx.clone();
        thread::spawn(move || {
            for row in (0..rows).filter(|r| r % num_procs == tid) {
                let mut vec_twi = vec![nodata; columns as usize];
                for col in 0..columns {
                    let z = dem.get_value(row, col);
                    if z != nodata {
                        let mut slope = match slope_type {
                            0_isize => get_gradient(&dem, row, col, z, cellsize), // local slope
                            1_isize => m_slope.get_value(row, col), // catchment slope
                            _ => panic!("Invalid 'slope_type' parameter"),
                        };

                        let slope2 = slope + slope_offset_rad;
                        slope = if slope2 > slope_min_rad { slope2.atan() } else { slope_min_rad.atan() };
        
                        let area = match area_type {
                            0_isize => m_amod.get_value(row, col), // total catchment area
                            1_isize => m_amod.get_value(row, col).sqrt(), // square root of catchment area
                            2_isize => m_amod.get_value(row, col) / cellsize as f32, // specific catchment area
                            _ => panic!("Invalid 'area_type' parameter"),
                        };
        
                        vec_twi[col as usize] = (area / slope).ln() as f64;
                    }
                }
                tx.send((row, vec_twi)).unwrap();
            }
        });
    }

    for _ in 0..rows {
        let (row, vec_twi) = rx.recv().expect("Error receiving data from thread.");
        twi.set_row_data(row, vec_twi);
    }
}
