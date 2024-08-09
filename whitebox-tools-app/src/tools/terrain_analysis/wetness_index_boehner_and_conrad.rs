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


/// This tool can be used to calculate the topographic wetness index, commonly used in the TOPMODEL rainfall-runoff framework.
/// The index describes the propensity for a site to be saturated to the surface given its contributing area and local slope
/// characteristics. It is calculated as:
///
/// > WI = Ln(As / tan(Slope))
///
/// Where `As` is the specific catchment area (i.e. the upslope contributing area per unit contour length) estimated using one of
/// the available flow accumulation algorithms in the Hydrological Analysis toolbox. Notice that `As` must not be log-transformed
/// prior to being used; log-transformation of `As` is a common practice when visualizing the data. The slope image should be
/// measured in degrees and can be created from the base digital elevation model (DEM) using the `Slope` tool. Grid cells with a
/// slope of zero will be assigned **NoData** in the output image to compensate for the fact that division by zero is infinity.
/// These very flat sites likely coincide with the wettest parts of the landscape. The input images must have the same grid dimensions.
///
/// Grid cells possessing the NoData value in either of the input images are assigned NoData value in the output image. The output
/// raster is of the float data type and continuous data scale.
/// 
/// Derived from the C++ implementation by Olaf Conrad for SAGA GIS 
///
/// See Also
/// `Slope`, `D8FlowAccumulation`, `DInfFlowAccumulation`, `FD8FlowAccumulation`, `BreachDepressionsLeastCost`, `WetnessIndex`
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
        let mut output_file = String::new();

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


            if vec[0].to_lowercase() == "-dem" || vec[0].to_lowercase() == "--dem" || vec[0].to_lowercase() == "-i"|| vec[0].to_lowercase() == "--input" {
                if keyval {
                    dem_file = vec[1].to_string();
                } else {
                    dem_file = args[i + 1].to_string();
                }
            } else if vec[0].to_lowercase() == "-o" || vec[0].to_lowercase() == "--output" {
                if keyval {
                    output_file = vec[1].to_string();
                } else {
                    output_file = args[i + 1].to_string();
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

        if !output_file.contains(&sep) && !output_file.contains("/") {
            output_file = format!("{}{}", working_directory, output_file);
        }
        if !dem_file.contains(&sep) && !dem_file.contains("/") {
            dem_file = format!("{}{}", working_directory, dem_file);
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
        let cellsize = resx;
        let num_procs = num_cpus::get() as isize;

        // À mettre en paramètres pour l'utilisateur
        let suction_ini = 10_f32;
        let slope_weight = 1_f32;
        let area_type = 0_isize;    // 0 = "total catchment area" | 1 = "square root of catchment area" | 2 = "specific catchment area"
        let slope_type = 1_isize;   // 0 = "local slope" | 1 = "catchment slope"
        let slope_min = 0_f32;      // degrees
        let slope_offset = 0_f32;   // degrees
        let mfd_converge = 1.1_f64;
        

        // Make sure that the DEM has square pixels
        // À modifier éventuellement pour permettre des pixels rectangulaires, je dois
        // juste m'assurer de ne pas me tromper d'où mettre resx et resy
        if resx != resy {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "The input DEM must have square pixels.",
            ));
        }





        ///////////
        // Calcul des accumulations de flux selon MFD
        // DEBUT DE get_area()
        ///////////
        println!("Calculate initial slope and suction matrices...");
        let mut m_suction: Array2D<f32> = Array2D::new(rows, columns, 0f32, -1f32)?;
        let mut m_slope: Array2D<f32> = Array2D::new(rows, columns, 0f32, -1f32)?;
        
        
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
                            let t_param = suction_ini.powf(slope_weight * slope);
                            vec_suction[col as usize] = (1.0 / t_param).powf(t_param.exp());

                            // Ajout à l'index d'élévation
                            vec_cells.push((row, col, z));

                        }
                    }
                    tx.send((row, vec_slope, vec_suction, vec_cells)).unwrap();
                }
            });
        }


        let mut cells_ordered = Vec::<(isize, isize, f64)>::with_capacity((rows * columns) as usize);
        for _ in 0..rows {
            let (row, vec_slope, vec_suction, mut vec_cells) = rx.recv().expect("Error receiving data from thread.");
            m_slope.set_row_data(row, vec_slope);
            m_suction.set_row_data(row, vec_suction);
            cells_ordered.append(&mut vec_cells);
        }

        // In order to pop the values from highest to lowest, we need to sort them from lowest to highest.
        // To ensure constant ordering from run to run (due to multiprocessing), values are first sorted by row and column
        cells_ordered.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Equal));
        cells_ordered.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Equal));
        cells_ordered.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(Equal));





        // Création de la matrice d'accumulation et de la matrice de poids
        let mut m_area: Array2D<f32> = Array2D::new(rows, columns, 0f32, -1f32)?;
        let mut m_weight: Array2D<f32> = Array2D::new(rows, columns, 1f32, -1f32)?;



        // Calcul du MFD initial
        // N'est finalement pas parallélisable à cause de la modification
        // itérative de "m_area" et "m_slope"
        println!("Calculate initial MFD...");
        let mut row: isize;
        let mut col: isize;

        let mut slope: f32;
        let mut t_param: f32;
        let mut area: f32;

        let (mut d, mut z): (f64, f64);
        let dcol = [0, 1, 1, 1, 0, -1, -1, -1];
        let drow = [1, 1, 0, -1, -1, -1, 0, 1];
        let diag_cellsize = (2.0 * cellsize * cellsize).sqrt();
        let grid_lengths = [
                    cellsize,
                    diag_cellsize,
                    cellsize,
                    diag_cellsize,
                    cellsize,
                    diag_cellsize,
                    cellsize,
                    diag_cellsize,
                ];


        while let Some(cell) = cells_ordered.pop() {
            row = cell.0;
            col = cell.1;
            z = cell.2;

            // Ajustement initial de l'accumulation (non parallélisable à cause de la modification itérative de "m_area" plus loin)
            area = m_area.get_value(row, col) + m_weight.get_value(row, col);
            m_area.set_value(row, col, area);

            // Ajustement initial de la pente du catchment en fonction de l'accumulation
            slope = m_slope.get_value(row, col);
            m_slope.set_value(row, col, slope / area);


            // Ajustement final de l'accumulation et de la pente du catchment
            let mut dz = vec![0_f32; 8];
            let mut dz_sum = 0_f32;
            let mut row_n: isize;
            let mut col_n: isize;
            let mut z_n: f64;

            for ii in 0..8 {
                row_n = row + drow[ii];
                col_n = col + dcol[ii];
                z_n = dem.get_value(row_n, col_n);
                if z_n != nodata {
                    d = z - z_n;
                    if d > 0.0 {
                        dz[ii] = (d / grid_lengths[ii]).atan().powf(mfd_converge) as f32;
                        dz_sum += dz[ii];
                    }
                }
            }

            if dz_sum > 0.0 {
                for ii in 0..8 {
                    if dz[ii] > 0.0 {
                        row_n = row + drow[ii];
                        col_n = col + dcol[ii];
                        z_n = dem.get_value(row_n, col_n);
                        if z_n != nodata {
                            m_area.increment(row_n, col_n, area * dz[ii] / dz_sum);
                            m_slope.increment(row_n, col_n, slope * dz[ii] / dz_sum);
                        }
                    }
                }
            }
        }


        // Ajustement de l'accumulation en fonction de la taille de cellule
        // Pas pertinent à paralléliser, il risque d'y avoir utimement plus de visites
        // de cellules
        let cell_area = (cellsize * cellsize) as f32;
        for row in 0..rows {
            for col in 0..columns {
                z = dem.get_value(row, col);
                if z != nodata {
                    m_area.set_value(row, col, m_area.get_value(row, col) * cell_area);
                }
            }
        }

        let mut raster_m_slope = Raster::initialize_using_array2d("Y:/Developpement/Programmation/JFB/Saga_TWI/m_slope.sdat", &dem.configs, m_slope.duplicate());
        let _ = raster_m_slope.write();

        let mut raster_m_suction = Raster::initialize_using_array2d("Y:/Developpement/Programmation/JFB/Saga_TWI/m_suction.sdat", &dem.configs, m_suction.duplicate());
        let _ = raster_m_suction.write();

        let mut raster_m_area = Raster::initialize_using_array2d("Y:/Developpement/Programmation/JFB/Saga_TWI/m_area.sdat", &dem.configs, m_area.duplicate());
        let _ = raster_m_area.write();

        // FIN DE get_area()




        // Création du masque de non traitement et ajustement des valeurs d'accumulation en
        // fonction de la couche d'accumulation utilisateur
        let mut m_mask: Array2D<i8> = Array2D::new(rows, columns, 0i8, -1i8)?; // Pour l'instant je traite partout (aucun masque)


        // Calcul itératif de l'accumulation de flux modifiée
        println!("Modify MFD...");
        let m_amod = get_modified(&m_area, m_suction, m_mask);

        let mut raster_m_amod = Raster::initialize_using_array2d("Y:/Developpement/Programmation/JFB/Saga_TWI/m_amod.sdat", &dem.configs, m_amod.duplicate());
        let _ = raster_m_amod.write();

        // Calcul classique du TWI
        println!("Calculate topographic wetness index...");
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
    // aspect = (-h/-g).atan()

    return slope;
}


fn get_modified<'a>(m_area_ini: &'a Array2D<f32>, m_suction: Array2D<f32>, m_mask: Array2D<i8>) -> Array2D<f32> {
    let rows = m_area_ini.rows as isize;
    let columns = m_area_ini.columns as isize;

    let mut masked: i8;  // Est-ce que je pourrais le remplir pour éviter des calculs pour les cellules nodata?
    let (mut area, mut area_mod): (f32, f32);

    let mut m_amod = m_area_ini.duplicate();
    let mut m_area = m_area_ini.duplicate();

    let mut nb_changes = 1usize;
    let mut iteration = 0usize;

    while nb_changes > 0 {
        iteration += 1;
        nb_changes = 0;

        // Boucle parallélisable... mais pas certain parce que "area_mod" nécessite
        // les valeurs voisine dans get_local_maximum() et que m_area est à la fois
        // lu et modifier. À valider dans le code C++. Ultimement, ce n'est peut-être pas
        // grave s'il y a convergence. À tester comme il faut.
        // Autre chose à tester: intégrer directement ici la fonction get_local_maximum
        // pour voir s'il y a un impact sur la performance (éviterait une déclaration répétitive
        // de variables, mais le compileur voit peut-être les choses autrement)
        for row in 0..rows {
            for col in 0..columns {
                // Un masque permet d'éviter à l'algorithme de perdre son temps dans une zone certaine d'eau (1 == eau, 0 == terre)
                // Pas nécessairement la meilleure approche, je devrais peut-être plutôt modifier "m_area_ini" en entrée avec le masque
                // pour y inscrire une valeur élevée d'accumulation sous le masque
                masked = m_mask.get_value(row, col);
                if masked == 0 {
                    area_mod = m_suction.get_value(row, col) * get_local_maximum(&m_area, row, col);
                    area = m_area.get_value(row, col);
                    if area_mod > area {
                        nb_changes += 1;
                        m_area.set_value(row, col, area_mod);
                    }

                }
            }
        }

        if nb_changes > 0 {
            nb_changes = 0;

            // Boucle parallélisable
            for row in 0..rows {
                for col in 0..columns {
                    area = m_area.get_value(row, col);
                    if area != m_amod.get_value(row, col) {
                        nb_changes += 1;
                        m_amod.set_value(row, col, area);
                    }
                }
            }
        }


        println!("pass {} ({} > 0)", iteration, nb_changes);
    }

    //let mut raster_m_areamod = Raster::initialize_using_array2d("Y:/Developpement/Programmation/JFB/Saga_TWI/m_amod_hatif.sdat", &dem.configs, m_amod.duplicate());
    //let _ = raster_m_areamod.write();

    println!("\npost-processing...");

    let (mut row_n, mut col_n): (isize, isize);
    let nodata_area = 0_f32; // Techniquement, ce n'est pas du nodata, mais ça y correspond. Formater différemment éventuellement.
    let m_area_nodata = -1_f32;

    // Boucle parallélisable car "m_amod" est jamais lu, toujours jsute modifié
    for row in 0..rows {
        for col in 0..columns {
            if m_area_ini.get_value(row, col) != nodata_area {
                let mut area_modified = false;
                let mut n = 0isize;
                let mut z = 0f32;

                for drow in -1..2 {
                    row_n = row + drow;
                    for dcol in -1..2 {
                        col_n = col + dcol;
                        if m_area_ini.get_value(row_n, col_n) != nodata_area {
                            area = m_area.get_value(row_n, col_n);
                            if area > m_area_ini.get_value(row_n, col_n) {
                                area_modified = true;
                            }
                            n += 1;
                            z += area;
                        }
                    }
                }
                if area_modified {
                    m_amod.set_value(row, col, z / n as f32)
                } else {
                    m_amod.set_value(row, col, m_area.get_value(row, col))
                }
            } else {
                m_amod.set_value(row, col, m_area_nodata)
            }

        }
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
    let cellsize = dem.configs.resolution_x;

    let (mut slope, mut slope2): (f32, f32);
    let mut area: f32;
    let mut z: f64;

    let slope_min_rad = slope_min * f32::consts::PI / 180.0;
    let slope_offset_rad = slope_offset * f32::consts::PI / 180.0;

    // Étape se parallélisant bien
    for row in 0..rows {
        for col in 0..columns {
            z = dem.get_value(row, col);
            if z != nodata {
                if slope_type == 1 {
                    // slope_type is "catchment slope"
                    slope = m_slope.get_value(row, col);
                } else {
                    // slope_type is "local slope"
                    slope = get_gradient(&dem, row, col, z, cellsize);
                }

                slope2 = slope + slope_offset_rad;
                slope = if slope2 > slope_min_rad { slope2.atan() } else { slope_min_rad.atan() };
                area = m_amod.get_value(row, col); //equivalent to area_type == "total catchment area"

                if area_type == 1isize {
                    // area_type is "square root of catchment area"
                    area = area.sqrt(); 
                } else if area_type == 2isize {
                    // area_type is "specific catchment area"
                    area = area / cellsize as f32;
                }

                twi.set_value(row, col, (area / slope).ln() as f64)
            }
        }
    }
}
