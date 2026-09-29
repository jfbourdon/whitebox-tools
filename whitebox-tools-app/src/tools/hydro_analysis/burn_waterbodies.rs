/*
This tool is part of the WhiteboxTools geospatial analysis library.
Authors: Dr. John Lindsay
Created: 17/04/2018
Last Modified: 18/10/2019
License: MIT
*/

use whitebox_common::algorithms::point_in_poly;
use whitebox_raster::*;
use whitebox_common::structures::{Array2D, BoundingBox, Point2D};
use crate::tools::*;
use whitebox_vector::{FieldData, ShapeType, Shapefile};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, VecDeque};
use std::env;
use std::f64;
use std::io::{Error, ErrorKind};
use std::path;

pub struct BurnWaterbodies {
    name: String,
    description: String,
    toolbox: String,
    parameters: Vec<ToolParameter>,
    example_usage: String,
}

impl BurnWaterbodies {
    /// public constructor
    pub fn new() -> BurnWaterbodies {
        let name = "BurnWaterbodies".to_string();
        let toolbox = "Hydrological Analysis".to_string();
        let description = "Burn waterbody polygons into a dem while enforcing flowpath.".to_string();

        let mut parameters = vec![];
        parameters.push(ToolParameter {
            name: "Input Waterbody Vector Polygon File".to_owned(),
            flags: vec!["--waterbodies".to_owned()],
            description: "Input waterbody vector polygons file.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Vector(
                VectorGeometryType::Polygon,
            )),
            default_value: None,
            optional: false,
        });

        parameters.push(ToolParameter {
            name: "Input Outlet Vector Point File".to_owned(),
            flags: vec!["--outlets".to_owned()],
            description: "Input outlet vector points file.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Vector(
                VectorGeometryType::Point,
            )),
            default_value: None,
            optional: false,
        });

        parameters.push(ToolParameter {
            name: "Input DEM File".to_owned(),
            flags: vec!["--dem".to_owned()],
            description: "Input raster DEM file.".to_owned(),
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
        let usage = format!(">>.*{0} -r={1} -v --wd=\"*path*to*data*\" -i=lakes.shp --field=ELEV -o=output.tif --nodata --cell_size=10.0
        >>.*{0} -r={1} -v --wd=\"*path*to*data*\" -i=lakes.shp --field=ELEV -o=output.tif --base=existing_raster.tif", short_exe, name).replace("*", &sep);

        BurnWaterbodies {
            name: name,
            description: description,
            toolbox: toolbox,
            parameters: parameters,
            example_usage: usage,
        }
    }
}

impl WhiteboxTool for BurnWaterbodies {
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
        let mut input_waterbodies = String::new();
        let mut input_outlets = String::new();
        let mut input_dem = String::new();
        let mut flat_increment = f64::NAN;
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
            let flag_val = vec[0].to_lowercase().replace("--", "-");
            if flag_val == "-i" || flag_val == "-waterbodies" {
                input_waterbodies = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            } else if flag_val == "-outlets" {
                input_outlets = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            } else if flag_val == "-dem" {
                input_dem = if keyval {
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

        if !input_waterbodies.contains(&sep) && !input_waterbodies.contains("/") {
            input_waterbodies = format!("{}{}", working_directory, input_waterbodies);
        }
        if !input_outlets.contains(&sep) && !input_outlets.contains("/") {
            input_outlets = format!("{}{}", working_directory, input_outlets);
        }
        if !input_dem.contains(&sep) && !input_dem.contains("/") {
            input_dem = format!("{}{}", working_directory, input_dem);
        }
        if !output_file.contains(&sep) && !output_file.contains("/") {
            output_file = format!("{}{}", working_directory, output_file);
        }


        if verbose {
            println!("Reading data...")
        };

        let start = Instant::now();


        // READING WATERBODIES VECTOR DATA
        let waterbodies_data = Shapefile::read(&input_waterbodies)?;

        // make sure the input vector file is of polygon type
        if waterbodies_data.header.shape_type.base_shape_type() != ShapeType::Polygon {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "The waterbodies vector data must be of polygon base shape type.",
            ));
        }



        // READING OUTLETS VECTOR DATA
        let outlets_data = Shapefile::read(&input_outlets)?;

        // make sure the input vector file is of polygon type
        if outlets_data.header.shape_type.base_shape_type() != ShapeType::Point {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "The oulets vector data must be of point base shape type.",
            ));
        }




        // Create the output raster. The process of doing this will
        // depend on whether a cell size or a base raster were specified.
        // If both are specified, the base raster takes priority.
        let mut dem = Raster::new(&input_dem, "r")?;
        //let mut output = Raster::initialize_using_file(&output_file, &dem);

        let rows = dem.configs.rows as isize;
        let columns = dem.configs.columns as isize;
        let raster_bb = dem.get_bounding_box();


        //////
        // Rasterisation initiale des polygones
        //////
        let mut start_point_in_part: usize;
        let mut end_point_in_part: usize;
        let mut col: isize;
        let mut row: isize;
        let mut output_something = false;
        let (mut x, mut y): (f64, f64);
        let (mut starting_row, mut ending_row, mut starting_col, mut ending_col): (
            isize,
            isize,
            isize,
            isize,
        );


        //// EMPÊCHER LA RASTERISATION DANS LE NODATA DU DEM !!!!
        /// À AJOUTER!!!  JE CROIS QUE ÇA VA RÉGLER MON PROBLÈME AVEC LES VALEURS TROP FAIBLES DANS LES SURFACES QUI DÉPASSENT

        let mut rasterized: Array2D<i8> = Array2D::new(rows, columns, -1i8, -1i8)?;
        let mut holes: Array2D<i32> = Array2D::new(rows, columns, -1i32, -1i32)?;
        let mut record_i32: i32;
        let num_records = waterbodies_data.num_records;
        for record_num in 0..waterbodies_data.num_records {
            let record = waterbodies_data.get_record(record_num);
            record_i32 = (record_num + 1) as i32;
            let rec_bb = BoundingBox::new(record.x_min, record.x_max, record.y_min, record.y_max);
            if rec_bb.overlaps(raster_bb) {
                // first find the holes
                for part in 0..record.num_parts as usize {
                    if record.is_hole(part as i32) {
                        start_point_in_part = record.parts[part] as usize;
                        end_point_in_part = if part < record.num_parts as usize - 1 {
                            record.parts[part + 1] as usize - 1
                        } else {
                            record.num_points as usize - 1
                        };

                        // First, figure out the minimum and maximum row and column for the polygon part
                        starting_row = rows;
                        ending_row = 0;
                        starting_col = columns;
                        ending_col = 0;
                        for p in start_point_in_part..end_point_in_part + 1 {
                            row = dem.get_row_from_y(record.points[p].y);
                            col = dem.get_column_from_x(record.points[p].x);
                            if row < starting_row {
                                starting_row = row;
                            }
                            if row > ending_row {
                                ending_row = row;
                            }
                            if col < starting_col {
                                starting_col = col;
                            }
                            if col > ending_col {
                                ending_col = col;
                            }
                        }

                        if starting_row < 0 {
                            starting_row = 0;
                        }
                        if ending_row < 0 {
                            ending_row = 0;
                        }
                        if starting_row >= rows {
                            starting_row = rows - 1;
                        }
                        if ending_row >= rows {
                            ending_row = rows - 1;
                        }

                        if starting_col < 0 {
                            starting_col = 0;
                        }
                        if ending_col < 0 {
                            ending_col = 0;
                        }
                        if starting_col >= columns {
                            starting_col = columns - 1;
                        }
                        if ending_col >= columns {
                            ending_col = columns - 1;
                        }

                        for r in starting_row..=ending_row {
                            y = dem.get_y_from_row(r);
                            for c in starting_col..=ending_col {
                                if dem.get_value(r, c) != dem.configs.nodata {
                                    x = dem.get_x_from_column(c);
                                    if point_in_poly(
                                        &Point2D { x: x, y: y },
                                        &record.points[start_point_in_part..=end_point_in_part],
                                    ) {
                                        holes.set_value(r, c, record_i32);
                                    }
                                }
                            }
                        }
                    }
                }

                for part in 0..record.num_parts as usize {
                    if !record.is_hole(part as i32) {
                        start_point_in_part = record.parts[part] as usize;
                        end_point_in_part = if part < record.num_parts as usize - 1 {
                            record.parts[part + 1] as usize - 1
                        } else {
                            record.num_points as usize - 1
                        };

                        // First, figure out the minimum and maximum row and column for the polygon part
                        starting_row = rows;
                        ending_row = 0;
                        starting_col = columns;
                        ending_col = 0;
                        for p in start_point_in_part..=end_point_in_part {
                            row = dem.get_row_from_y(record.points[p].y);
                            col = dem.get_column_from_x(record.points[p].x);
                            if row < starting_row {
                                starting_row = row;
                            }
                            if row > ending_row {
                                ending_row = row;
                            }
                            if col < starting_col {
                                starting_col = col;
                            }
                            if col > ending_col {
                                ending_col = col;
                            }
                        }

                        if starting_row < 0 {
                            starting_row = 0;
                        }
                        if ending_row < 0 {
                            ending_row = 0;
                        }
                        if starting_row >= rows {
                            starting_row = rows - 1;
                        }
                        if ending_row >= rows {
                            ending_row = rows - 1;
                        }

                        if starting_col < 0 {
                            starting_col = 0;
                        }
                        if ending_col < 0 {
                            ending_col = 0;
                        }
                        if starting_col >= columns {
                            starting_col = columns - 1;
                        }
                        if ending_col >= columns {
                            ending_col = columns - 1;
                        }

                        for r in starting_row..=ending_row {
                            y = dem.get_y_from_row(r);
                            for c in starting_col..=ending_col {
                                if dem.get_value(r, c) != dem.configs.nodata {
                                    x = dem.get_x_from_column(c);
                                    if point_in_poly(
                                        &Point2D { x: x, y: y },
                                        &record.points[start_point_in_part..end_point_in_part + 1],
                                    ) {
                                        if holes.get_value(r, c) != record_i32 {
                                            rasterized.set_value(r, c, 1i8);
                                            output_something = true;
                                        }
                                    }
                                }
                            }
                        }


                    // find each intersection of the polyline with a row
                    for r in starting_row..=ending_row {
                        let row_y_coord = dem.get_y_from_row(r);
                        // find the x-coordinates of each of the line segments
                        // that intersect this row's y coordinate
                        for i in start_point_in_part..end_point_in_part {
                            if is_between(row_y_coord, record.points[i].y, record.points[i + 1].y) {
                                let y1 = record.points[i].y;
                                let y2 = record.points[i + 1].y;
                                if y2 != y1 {
                                    let x1 = record.points[i].x;
                                    let x2 = record.points[i + 1].x;

                                    // calculate the intersection point
                                    let x_prime = x1 + (row_y_coord - y1) / (y2 - y1) * (x2 - x1);
                                    let c = dem.get_column_from_x(x_prime);
                                    if dem.get_value(r, c) != dem.configs.nodata {
                                        rasterized.set_value(r, c, 1i8);
                                        output_something = true;
                                    }
                                }
                            }
                        }
                    }

                    // find each intersection of the polyline with a column
                    for c in starting_col..=ending_col {
                        let col_x_coord = dem.get_x_from_column(c);
                        for i in start_point_in_part..end_point_in_part {
                            if is_between(col_x_coord, record.points[i].x, record.points[i + 1].x) {
                                let x1 = record.points[i].x;
                                let x2 = record.points[i + 1].x;
                                if x1 != x2 {
                                    let y1 = record.points[i].y;
                                    let y2 = record.points[i + 1].y;

                                    // calculate the intersection point
                                    let y_prime = y1 + (col_x_coord - x1) / (x2 - x1) * (y2 - y1);
                                    let r = dem.get_row_from_y(y_prime);
                                    if dem.get_value(r, c) != dem.configs.nodata {
                                        rasterized.set_value(r, c, 1i8);
                                        output_something = true;
                                    }
                                }
                            }
                        }
                    }




                    }
                }
            }
            if verbose {
                progress = (100.0_f64 * (record_num + 1) as f64 / num_records as f64) as usize;
                //if progress != old_progress {
                println!(
                    "Rasterizing {} of {}: {}%",
                    record_num + 1,
                    num_records,
                    progress
                );
                //    old_progress = progress;
                //}
            }
        }


        ////
        /// Nombre de cellules uniques visitées
        let mut nb_visited_cells = 0_usize;
        for row in 0..rows {
            for col in 0..columns {
                let z = rasterized.get_value(row, col);
                if z != rasterized.nodata {
                    nb_visited_cells += 1_usize
                }
            }
        }



        //////
        /// Début de l'opération de region-growing à partir d'un point d'exutoire
        //////
        let mut visited: Array2D<u8> = Array2D::new(rows, columns, 0u8, 0u8)?;
        let dx = [1, 1, 1, 0, -1, -1, -1, 0];
        let dy = [-1, 0, 1, 1, 1, 0, -1, -1];


        let num_records = outlets_data.num_records;
        for record_num in 0..outlets_data.num_records {
            let record = outlets_data.get_record(record_num);


            println!("record_num: {}", record_num);



            // Ajout du point comme point d'initiation de la recherche
            // Les multipoints permettent de partir en même temps de deux endroits différents
            // sur une même surface hydrographique.
            // Pas besoin d'un maxheap en fait, un Vec classé à la toute fin comme "undefined_flow_cells"
            // dans breach_depressions_least_cost ferait la même chose (possiblement à coût moindre en plus)
            // Pour déterminer la taille du max heap, mettre un compteur dans la phase de rasterisation
            // et sauvegarder le nombre total de cellules modifiées. En fait, mieux de repasser à travers
            // "rasterized" et counter le nombre différent de nodata car sinon j'aurais des superpositions.
            let mut next_queue: VecDeque<(isize, isize, usize)> = VecDeque::with_capacity(nb_visited_cells);
            //let mut maxheap  = BinaryHeap::with_capacity(nb_visited_cells);
            let mut maxheap: BinaryHeap<(GridCell)> = BinaryHeap::with_capacity(nb_visited_cells);

            for no_pt in 0..record.num_points as usize {
                let row_ini = dem.get_row_from_y(record.points[no_pt].y);
                let col_ini = dem.get_column_from_x(record.points[no_pt].x);
                next_queue.push_back((row_ini, col_ini, 0_usize));
                visited.set_value(row_ini, col_ini, 1u8);
            }


            // region growing operation to find all attached valid cells
            while !next_queue.is_empty() {
                let cell = next_queue.pop_front().unwrap();
                maxheap.push(GridCell {
                    row: cell.0,
                    column: cell.1,
                    priority: cell.2,
                });
                let next_priority = cell.2 + 1_usize;
                for n in 0..8 {
                    let row = cell.0 + dy[n];
                    let col = cell.1 + dx[n];
                    if rasterized.get_value(row, col) != rasterized.nodata &&  // Voir à combiner "rasterized" et "visited"
                    visited.get_value(row, col) == visited.nodata {
                        next_queue.push_back((row, col, next_priority));
                        visited.set_value(row, col, 1u8);
                    }
                }
            }


            // Poper "maxheap" en partant du point le plus éloigné de l'exutoire (en terme d'itération)
            println!("Longueur du maxheap: {}", maxheap.len());
            while !maxheap.is_empty() {
                let cell = maxheap.pop().expect("Error during pop operation.");
                let row = cell.row;
                let col = cell.column;
                let mut z = dem.get_value(row, col);

                // Recherche de la valeur d'élévation la plus faible parmi les voisins
                for n in 0..8 {
                    let row2 = row + dy[n];
                    let col2 = col + dx[n];
                    let zn = dem.get_value(row2, col2);
                    if zn != dem.configs.nodata && zn < z {
                        z = zn;
                    }
                }
                dem.set_value(row, col, z - flat_increment);
            }
        }



        /* 
        let s : &str = "F:\\06aj_026\\visited.tif";
        let mut output_raster = Raster::initialize_using_file(&s, &dem);
        output_raster.configs.data_type = DataType::U32;
        for row in 0..rows {
            for col in 0..columns {
                output_raster.set_value(row, col, visited.get_value(row, col) as f64);
            }
        }
        let _ = match output_raster.write() {
            Ok(_) => {
                if verbose {
                    println!("Output output_raster VISITED")
                }
            }
            Err(e) => return Err(e),
        };
        */



        // Redirect DEM raster to the output file
        // Force F64 to ensure correct increment values
        dem.file_name = output_file.to_owned();
        dem.file_mode = "w".to_string();
        dem.raster_type = get_raster_type_from_file(output_file.to_owned(), "w".to_string());
        dem.configs.data_type = DataType::F64;



        let elapsed_time = get_formatted_elapsed_time(start);
        dem.add_metadata_entry(format!(
            "Created by whitebox_tools\' {} tool",
            self.get_tool_name()
        ));
        dem.add_metadata_entry(format!("Elapsed Time (excluding I/O): {}", elapsed_time));

        if verbose {
            println!("Saving data...")
        };
        let _ = match dem.write() {
            Ok(_) => {
                if verbose {
                    println!("Output file written")
                }
            }
            Err(e) => return Err(e),
        };

        if !output_something && verbose {
            println!("Warning: No polygons were output to the raster.");
        }

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
    priority: usize,
}

impl Eq for GridCell {}

impl PartialOrd for GridCell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.priority.partial_cmp(&other.priority)
    }
}

impl Ord for GridCell {
    fn cmp(&self, other: &GridCell) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}


fn is_between(val: f64, threshold1: f64, threshold2: f64) -> bool {
    if val == threshold1 || val == threshold2 {
        return true;
    }
    if threshold2 > threshold1 {
        return val > threshold1 && val < threshold2;
    }
    val > threshold2 && val < threshold1
}