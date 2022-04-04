/*
This tool is part of the WhiteboxTools geospatial analysis library.
Authors: Jean-François Bourdon
Created: 04/04/2022
Last Modified: 04/04/2022
License: MIT
*/

// IMPORT généraux
use whitebox_raster::*;
use crate::tools::*;
use std::env;
use std::f64;
use std::io::{Error, ErrorKind};
use std::path;


// IMPORT vector_lines_to_raster.rs
use whitebox_common::structures::BoundingBox;
use whitebox_vector::{FieldData, ShapeType, Shapefile};


// IMPORT cost_distance.rs
use whitebox_common::structures::Array2D;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::i32;




/// This tool decrements (lowers) the elevations of pixels within an input digital elevation model (DEM) (`--dem`)
/// along an input vector stream network (`--streams`). The optional flat height increment value (`--flat_increment`)
/// can be specified to control the pixel lowering rate. Notice that **if the `--flat_increment` parameter is not specified,
/// the small number used to ensure flow across flats will be calculated automatically, which should be preferred in
/// most applications** of the tool. Stream vectors can be processed in a specific order by providing to `--order_by`
/// parameter the name of a numeric (integer) field from the attributes table. Processing is done in increasing order and
/// follow the numerisation direction of each stream vector to decrement elevations (from upstream to downstream).
///
/// All stream vectors must only cover valid cell (no NoData) from the DEM otherwise the tool is halted.
///
/// # See Also
/// `BurnStreamsAtRoads`, `RasterStreamsToVector`, `RasterizeStreams`
pub struct BurnStreams {
    name: String,
    description: String,
    toolbox: String,
    parameters: Vec<ToolParameter>,
    example_usage: String,
}

impl BurnStreams {
    /// public constructor
    pub fn new() -> BurnStreams {
        let name = "BurnStreams".to_string();
        let toolbox = "Hydrological Analysis".to_string();
        let description = "Burn-in streams into an elevation raster.".to_string();

        let mut parameters = vec![];
        parameters.push(ToolParameter{
            name: "Input DEM File".to_owned(),
            flags: vec!["--dem".to_owned()],
            description: "Input raster digital elevation model (DEM) file.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Raster),
            default_value: None,
            optional: false
        });

        parameters.push(ToolParameter {
            name: "Input Vector Streams File".to_owned(),
            flags: vec!["--streams".to_owned()],
            description: "Input vector streams file.".to_owned(),
            parameter_type: ParameterType::ExistingFile(ParameterFileType::Vector(
                VectorGeometryType::Line,
            )),
            default_value: None,
            optional: false,
        });

        parameters.push(ToolParameter {
            name: "Field Name".to_owned(),
            flags: vec!["--order_by".to_owned()],
            description: "Input field name in attributes table by which prioritise stream burning.".to_owned(),
            parameter_type: ParameterType::VectorAttributeField(
                AttributeType::Number,
                "--streams".to_string(),
            ),
            default_value: Some("FID".to_owned()),
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
        let usage = format!(">>.*{0} -r={1} -v --wd=\"*path*to*data*\" --dem=raster.tif --streams=lines.shp -o=output.tif
        >>.*{0} -r={1} -v --wd=\"*path*to*data*\" --dem=raster.tif --streams=lines.shp -o=output.tif --order_by=PRIORITY --flat_increment=0.01", short_exe, name).replace("*", &sep);

        BurnStreams {
            name: name,
            description: description,
            toolbox: toolbox,
            parameters: parameters,
            example_usage: usage,
        }
    }
}

impl WhiteboxTool for BurnStreams {
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
        let mut dem_file = String::new();
        let mut streams_file = String::new();
        let mut field_name = String::from("FID");
        let mut output_file = String::new();
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
            if flag_val == "-dem" {
                dem_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            } else if flag_val == "-streams" {
                streams_file = if keyval {
                    vec[1].to_string()
                } else {
                    args[i + 1].to_string()
                };
            } else if flag_val == "-order_by" {
                field_name = if keyval {
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

        if !streams_file.contains(&sep) && !streams_file.contains("/") {
            streams_file = format!("{}{}", working_directory, streams_file);
        }
        if !output_file.contains(&sep) && !output_file.contains("/") {
            output_file = format!("{}{}", working_directory, output_file);
        }

        if verbose {
            println!("Reading data...")
        };
        let streams = Shapefile::read(&streams_file).expect("Error reading input Shapefile.");

        let start = Instant::now();

        // make sure the input vector file is of polyline type (eventually should explicitly
        // look for LINESTRING only and not MULTILINESTRING)
        if streams.header.shape_type.base_shape_type() != ShapeType::PolyLine
        {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "The input vector data must either be of polyline base shape type.",
            ));
        }

        // What is the index of the field to be analyzed?
        let field_index = match streams.attributes.get_field_num(&field_name) {
            Some(i) => i,
            None => {
                // Field not found use FID
                if verbose {
                    println!("Warning: Attribute not found in table. FID will be used instead.");
                }
                field_name = "FID".to_string();
                0
            }
        };

        // Is the field numeric?
        if !streams.attributes.is_field_numeric(field_index) {
            // Warn user of non-numeric
            if verbose {
                println!("Warning: Non-numeric attributes cannot be rasterized. FID will be used instead.");
            }
            field_name = "FID".to_string(); // Can't use non-numeric field; use FID instead.
        }

        // Load reference raster
        let mut dem = Raster::new(&dem_file, "r")?;
        
        let background_val = dem.configs.nodata;
        let rows = dem.configs.rows as isize;
        let columns = dem.configs.columns as isize;


        // Transformation du Array2D en Raster simplement pour pouvoir avoir accès à
        // la famille de méthodes get_row_from_y. Je pourrais probablement m'en passer en
        // utilisant des offsets à partir du raster dem
        let rasterization_file = format!("{}{}", working_directory, "rasterization_file.flt");
        let mut rasterization = Raster::initialize_using_config(&rasterization_file, &dem.configs);
        rasterization.reinitialize_values(background_val);
        //let mut rasterization: Array2D<f64> = Array2D::new(rows, columns, 0.0f64, -1.0f64)?;

        //let mut rasterization = Raster::initialize_using_array2d(&rasterization_file, &dem.configs, rasterization_array);
        //rasterization.configs.data_type = DataType::I8;
        //rasterization.configs.nodata = 0i8;

        
        let mut attribute_data_ori = vec![0u64; streams.num_records];
        // get the attribute data
        for record_num in 0..streams.num_records {
            if field_name != "FID" {
                match streams.attributes.get_value(record_num, &field_name) {
                    FieldData::Int(val) => {
                        attribute_data_ori[record_num] = val as u64;
                    }
                    _ => {
                        // do nothing; likely due to null value for record.
                    }
                }
            } else {
                attribute_data_ori[record_num] = (record_num + 1) as u64;
            }

            if verbose {
                progress =
                    (100.0_f64 * record_num as f64 / (streams.num_records - 1) as f64) as usize;
                if progress != old_progress {
                    println!("Reading attributes: {}%", progress);
                    old_progress = progress;
                }
            }
        }


        // CLASSEMENT DES VECTEURS SELON L'ORDRE DE PRIORITÉ
        let idx: Vec<usize> = (0..streams.num_records).collect();
        let attribute_data_iter = attribute_data_ori.into_iter();

        let mut tuple_attribute_data: Vec<(usize, u64)> = idx.into_iter().zip(attribute_data_iter).collect();

        //https://stackoverflow.com/questions/40091161/sorting-a-vector-of-tuples-needs-a-reference-for-the-second-value
        tuple_attribute_data.sort_by_key(|k| k.1);

        // JE DOIS EXTRAIRE ENSUITE LES IDX QUI SONT À LA POSITION 0 DE CHAQUE TUPLE
        // C'est dégeu mais ça marche
        let mut record_num_sorted = vec![0; streams.num_records];
        for record_num in 0..streams.num_records {
            record_num_sorted[record_num] = tuple_attribute_data[record_num].0;
        }
        

        let raster_bb = BoundingBox::new(
            dem.configs.west,
            dem.configs.east,
            dem.configs.south,
            dem.configs.north,
        );
        let mut bb = BoundingBox {
            ..Default::default()
        };
        let (mut top_row, mut bottom_row, mut left_col, mut right_col): (
            isize,
            isize,
            isize,
            isize,
        );
        let mut row_y_coord: f64;
        let mut col_x_coord: f64;
        let (mut x1, mut x2, mut y1, mut y2): (f64, f64, f64, f64);
        let (mut x_prime, mut y_prime): (f64, f64);
        let mut start_point_in_part: usize;
        let mut end_point_in_part: usize;
        let mut output_something = false;
        let num_records = streams.num_records;


        for record_num in record_num_sorted {
            let record = streams.get_record(record_num);


            // Ensure that all streams are entirely contained within the raster extent
            // I could use the extent in the Shapefile header but the values seem to only be updated
            // when a feature is added or removed (at least with editing with QGIS)
            let rec_bb = BoundingBox::new(record.x_min, record.x_max, record.y_min, record.y_max);
            if !rec_bb.within(raster_bb)
            {
                return Err(Error::new(
                    ErrorKind::InvalidInput,
                    "The input vector data must be contained inside the raster extent. Some lines are currently extending beyond.",
                ));
            }


            // Ensure that only singlepart polylines are used. This case should be catch right at loading
            // but Whitebox doesn't distinguish currently singlepart from multipart
            // (see lines 500 and up in whitebox-vector\src\shapefile\geometry.rs)
            if record.num_parts as usize > 1
            {
                return Err(Error::new(
                    ErrorKind::InvalidInput,
                    "The input vector data must be of polyline singlepart type, not polyline multipart.",
                ));
            }

            
            start_point_in_part = 0 as usize;
            end_point_in_part = record.num_points as usize - 1;

            bb.initialize_to_inf();
            for i in start_point_in_part..end_point_in_part + 1 {
                if record.points[i].x < bb.min_x {
                    bb.min_x = record.points[i].x;
                }
                if record.points[i].x > bb.max_x {
                    bb.max_x = record.points[i].x;
                }
                if record.points[i].y < bb.min_y {
                    bb.min_y = record.points[i].y;
                }
                if record.points[i].y > bb.max_y {
                    bb.max_y = record.points[i].y;
                }
            }
            top_row = dem.get_row_from_y(bb.max_y);
            bottom_row = dem.get_row_from_y(bb.min_y);
            left_col = dem.get_column_from_x(bb.min_x);
            right_col = dem.get_column_from_x(bb.max_x);

            if top_row < 0 {
                top_row = 0;
            }
            if bottom_row < 0 {
                bottom_row = 0;
            }
            if top_row >= rows {
                top_row = rows - 1;
            }
            if bottom_row >= rows {
                bottom_row = rows - 1;
            }

            if left_col < 0 {
                left_col = 0;
            }
            if right_col < 0 {
                right_col = 0;
            }
            if left_col >= columns {
                left_col = columns - 1;
            }
            if right_col >= columns {
                right_col = columns - 1;
            }
            
            // Trouve les points de départ et d'arrivée de la ligne sur la grille
            // Les points de départ et d'arrivée sont forcés dans le raster pour m'assurer
            // qu'ils soient toujours utilisés. Assure aussi un cohérance si des lignes
            // indépendantes sont connectées ensemble
            let row_start = dem.get_row_from_y(record.points[start_point_in_part].y);
            let col_start = dem.get_column_from_x(record.points[start_point_in_part].x);
            let row_end = dem.get_row_from_y(record.points[end_point_in_part].y);
            let col_end = dem.get_column_from_x(record.points[end_point_in_part].x);

            rasterization.set_value(row_start, col_start, 1.0f64);
            rasterization.set_value(row_end, col_end, 1.0f64);
            

            // find each intersection with a row.
            for row in top_row..bottom_row + 1 {
                row_y_coord = dem.get_y_from_row(row);
                // find the x-coordinates of each of the line segments
                // that intersect this row's y coordinate
                for i in start_point_in_part..end_point_in_part {
                    if is_between(row_y_coord, record.points[i].y, record.points[i + 1].y) {
                        y1 = record.points[i].y;
                        y2 = record.points[i + 1].y;
                        if y2 != y1 {
                            x1 = record.points[i].x;
                            x2 = record.points[i + 1].x;

                            // calculate the intersection point
                            x_prime = x1 + (row_y_coord - y1) / (y2 - y1) * (x2 - x1);
                            let col = dem.get_column_from_x(x_prime);

                            rasterization.set_value(row, col, 1.0f64);
                            output_something = true;
                        }

                    }
                }
            }

            // find each intersection with a column.
            for col in left_col..right_col + 1 {
                col_x_coord = dem.get_x_from_column(col);
                for i in start_point_in_part..end_point_in_part {
                    if is_between(col_x_coord, record.points[i].x, record.points[i + 1].x) {
                        x1 = record.points[i].x;
                        x2 = record.points[i + 1].x;
                        if x1 != x2 {
                            y1 = record.points[i].y;
                            y2 = record.points[i + 1].y;

                            // calculate the intersection point
                            y_prime = y1 + (col_x_coord - x1) / (x2 - x1) * (y2 - y1);
                            let row = dem.get_row_from_y(y_prime);

                            rasterization.set_value(row, col, 1.0f64);
                            output_something = true;
                        }
                    }
                }
            }

            //println!("Start: R{} C{}   End: R{} C{}", row_start, col_start, row_end, col_end);


            // UTILISATION DE COST DISTANCE POUR FAIRE LA MATRICE DE COUTS

            // Se trouve à être unique un raster indiquant la position de la cellule
            // de fin de la ligne. La cellule de fin est utilisée car la fonction effectue
            // un backlink. Les directions de flux trouvées indique donc le chemin vers la
            // cellule désignée faisant en sorte que cette façon de faire me permet d'avoir
            // le trajet du point départ vers le point d'arrivée et non l'inverse.
            let mut source: Array2D<i8> = Array2D::new(rows, columns, 0, 0)?;
            source.set_value(row_end, col_end, 1i8);

            let cost = rasterization.clone(); // se trouve à être l'output de la section précédente

            let num_cells = (rows * columns) as usize;
            let nodata = rasterization.configs.nodata;
            
            let accum_file = format!("{}{}", working_directory, "accum_file.flt");
            let mut output2 = Raster::initialize_using_file(&accum_file, &dem);
            output2.configs.data_type = DataType::F32;
            let background_val = (i32::max_value() - 1) as f64;
            output2.reinitialize_values(background_val);
            
            let backlink_file = format!("{}{}", working_directory, "backlink_file.flt");
            let mut backlink = Raster::initialize_using_file(&backlink_file, &dem);
    
            let mut minheap = BinaryHeap::with_capacity(num_cells);
    
            // Forte opimisation possible puisque je sais déjà quelle
            // cellule utiliser
            for row in 0..rows {
                for col in 0..columns {
                    if source.get_value(row, col) == 1i8 && cost.get_value(row, col) != nodata {
                        output2.set_value(row, col, 0.0);
                        backlink.set_value(row, col, 0.0);
                    } else if cost.get_value(row, col) == nodata {
                        output2.set_value(row, col, nodata);
                    }
                }
            }

            minheap.push(GridCell {
                row: row_end,
                column: col_end,
                priority: 0f64,
            });
    
            let mut new_cost: f64;
            let mut accum_val: f64;
            let (mut cost1, mut cost2): (f64, f64);
            let (mut row, mut col): (isize, isize);
            let (mut row_n, mut col_n): (isize, isize);
            let cell_size_x = dem.configs.resolution_x;
            let cell_size_y = dem.configs.resolution_y;
            let diag_cell_size = (cell_size_x * cell_size_x + cell_size_y * cell_size_y).sqrt();
            let dist = [
                diag_cell_size,
                cell_size_x,
                diag_cell_size,
                cell_size_y,
                diag_cell_size,
                cell_size_x,
                diag_cell_size,
                cell_size_y,
            ];
            let dx = [1, 1, 1, 0, -1, -1, -1, 0];
            let dy = [-1, 0, 1, 1, 1, 0, -1, -1];
            let backlink_dir = [16.0, 32.0, 64.0, 128.0, 1.0, 2.0, 4.0, 8.0];
            let mut solved: Array2D<i8> = Array2D::new(rows, columns, 0, -1)?;
            while !minheap.is_empty() {
                let cell = minheap.pop().expect("Error during pop operation.");
                row = cell.row;
                col = cell.column;
                if solved.get_value(row, col) == 0 {
                    solved.set_value(row, col, 1);
                    accum_val = output2.get_value(row, col);
                    cost1 = cost.get_value(row, col);
                    for n in 0..8 {
                        col_n = col + dx[n];
                        row_n = row + dy[n];
                        if output2.get_value(row_n, col_n) != nodata {
                            cost2 = cost.get_value(row_n, col_n);
                            new_cost = accum_val + (cost1 + cost2) / 2.0 * dist[n];
                            if new_cost < output2.get_value(row_n, col_n) {
                                if solved.get_value(row_n, col_n) == 0 {
                                    output2.set_value(row_n, col_n, new_cost);
                                    backlink.set_value(row_n, col_n, backlink_dir[n]);
                                    minheap.push(GridCell {
                                        row: row_n,
                                        column: col_n,
                                        priority: new_cost,
                                    });
                                }
                            }
                        }
                    }
                }
            }
            
            /*
            let _ = match output2.write() {
                Ok(_) => {
                    if verbose {
                        println!("Output accum_file file written")
                    }
                }
                Err(e) => return Err(e),
            };

            let _ = match backlink.write() {
                Ok(_) => {
                    if verbose {
                        println!("Output backlink file written")
                    }
                }
                Err(e) => return Err(e),
            };
            */

            // UTILISATION DU BACKLINK POUR FAIRE DÉCROITRE LES VALEURS D'ÉLÉVATION
            // REPREND DES SECTIONS DE breach_depressions.rs
            // part du début du segment et remonte le backlink jusqu'à 0
            let nodata = dem.configs.nodata;
            let resx = dem.configs.resolution_x;
            let resy = dem.configs.resolution_y;
            let diagres = (resx * resx + resy * resy).sqrt();


            let small_num = if !flat_increment.is_nan() || flat_increment == 0f64 {
                flat_increment
            } else {
                let elev_digits = (dem.configs.maximum as i64).to_string().len();
                let elev_multiplier = 10.0_f64.powi((6 - elev_digits) as i32);
                1.0_f64 / elev_multiplier as f64 * diagres.ceil()
            };


            let dx = [1, 1, 1, 0, -1, -1, -1, 0];
            let dy = [-1, 0, 1, 1, 1, 0, -1, -1];

            let mut d8pointer = backlink.get_value(row_start, col_start);
            let mut zin: f64;
            let mut z_target: f64;
            let mut zin_n: f64;
            let (mut row_n, mut col_n): (isize, isize);
            let mut dir: f64;

            row = row_start;
            col = col_start;

            while d8pointer != 0f64{
                //println!("--------\n- Pointer: {}   ROW: {}   COL: {}", d8pointer, row, col);
                zin = dem.get_value(row, col);
                z_target = zin;
                
                // Compare elevation of each neighbour to the central cell
                // and lower it to the bare minimum to clear the lower cells
                for n in 0..8 {
                    row_n = row + dy[n];
                    col_n = col + dx[n];
                    zin_n = dem.get_value(row_n, col_n);
                    //println!("-- ZIN_N: {}", zin_n);
                    if zin_n != nodata {
                        if zin_n <= z_target {
                            //println!("{} to {}", z_target, zin_n - small_num);
                            z_target = zin_n - small_num;
                        }
                    }
                }
                dem.set_value(row, col, z_target);
                //println!("-- ZIN: {}   Z_TARGET: {}\n", zin, z_target);

                // Get the next cell to compare
                // je pourrais me débarraser du log2 et de la condition si le backlink produisait
                // des valeurs de 0-1-2-3-4-5-6-7 plutôt que 0-2-4-8-16-32-64-128
                d8pointer = backlink.get_value(row, col);
                if d8pointer != 0f64 {
                    dir = d8pointer.log2();
                    row += dy[dir as usize];
                    col += dx[dir as usize];
                }
            }

            if verbose {
                progress = (100.0_f64 * (record_num + 1) as f64 / num_records as f64) as usize;
                if progress != old_progress {
                    println!(
                        "Rasterizing {} of {}: {}%",
                        record_num + 1,
                        num_records,
                        progress
                    );
                    old_progress = progress;
                }
            }
        }

        let elapsed_time = get_formatted_elapsed_time(start);


        let mut dem_updated = Raster::initialize_using_file(&output_file, &dem);
        // Réinscrit les valeurs du fichier de dem
        // Il y a probablement une manière plus intelligente de faire ça
        for row in 0..rows {
            let data = dem.get_row_data(row);
            dem_updated.set_row_data(row, data);
        }

        let _ = match dem_updated.write() {
            Ok(_) => {
                if verbose {
                    println!("Output burnt file written (dem_updated)")
                }
            }
            Err(e) => return Err(e),
        };


        if !output_something && verbose {
            println!("Warning: No polylines were output to the raster.");
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

/// FONCTIONS POUR SECTION vector_lines_to_raster.rs
fn is_between(val: f64, threshold1: f64, threshold2: f64) -> bool {
    if val == threshold1 || val == threshold2 {
        return true;
    }
    if threshold2 > threshold1 {
        return val > threshold1 && val < threshold2;
    }
    val > threshold2 && val < threshold1
}


/// FONCTIONS POUR SECTION cost_distance.rs
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
