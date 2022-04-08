/*
This tool is part of the WhiteboxTools geospatial analysis library.
Authors: Jean-François Bourdon
Created: 04/04/2022
Last Modified: 08/04/2022
License: MIT
*/

use whitebox_raster::*;
use whitebox_common::structures::BoundingBox;
use whitebox_common::structures::Array2D;
use crate::tools::*;
use whitebox_vector::{FieldData, ShapeType, Shapefile};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::env;
use std::i32;
use std::f64;
use std::io::{Error, ErrorKind};
use std::path;

/// This tool decrements (lowers) the elevation of pixels within an input digital elevation model (DEM) (`--dem`)
/// along an input vector stream network (`--streams`) following the numerisation direction of each stream vector
/// (from upstream to downstream). Streams can be processed in a specific order by providing to the `--order_by`
/// parameter the name of a numeric (integer) field from the attributes table (ascending order).
/// The optional flat height increment value (`--flat_increment`) can be specified to control
/// the pixel lowering rate. Notice that if the `--flat_increment` parameter isn't specified,
/// the output DEM will use 64-bit floating point values in order to make sure that the very small elevation increment
/// value determined will be accurately stored. Consequently, it may double the storage requirements as DEMs are often
/// stored with 32-bit precision. However, if a flat increment value is specified, the output DEM will keep the same data
/// type as the input assuming the user chose its value wisely.
///
/// Some restrictions to the stream vectors geometry are in place in order to prevent any guessing work. Streams must all
/// be single part and only cover valid cells (no NoData) from the DEM.
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
        let mut iter_progess: usize = 0;
        let mut old_progress: usize = 1;

        if !streams_file.contains(&sep) && !streams_file.contains("/") {
            streams_file = format!("{}{}", working_directory, streams_file);
        }
        if !output_file.contains(&sep) && !output_file.contains("/") {
            output_file = format!("{}{}", working_directory, output_file);
        }

        if verbose {
            println!("Reading stream data...")
        };
        let streams = Shapefile::read(&streams_file).expect("Error reading input Shapefile.");

        let start = Instant::now();

        // Make sure the input vector file is of polyline type
        // Should explicitly look for LINESTRING only and not MULTILINESTRING ideally
        if streams.header.shape_type.base_shape_type() != ShapeType::PolyLine
        {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "The input vector data must be of polyline base shape type.",
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
                println!("Warning: Non-numeric attributes cannot be used for ordering. FID will be used instead.");
            }
            field_name = "FID".to_string(); // Can't use non-numeric field; use FID instead.
        }


        // Get the attribute data
        let mut attribute_data = vec![0u64; streams.num_records];
        for record_num in 0..streams.num_records {
            if field_name != "FID" {
                match streams.attributes.get_value(record_num, &field_name) {
                    FieldData::Int(val) => {
                        attribute_data[record_num] = val as u64;
                    }
                    _ => {
                        // do nothing; likely due to null value for record.
                    }
                }
            } else {
                attribute_data[record_num] = (record_num + 1) as u64;
            }
        }


        // Sort stream vectors based on the priority field
        // My way of doing it is really ugly but it works
        // There is certainly a more elegent way of doing this
        let idx: Vec<usize> = (0..streams.num_records).collect();
        let attribute_data_iter = attribute_data.into_iter();
        let mut tuple_attribute_data: Vec<(usize, u64)> = idx.into_iter().zip(attribute_data_iter).collect();

        tuple_attribute_data.sort_by_key(|k| k.1);
        let mut record_num_sorted = vec![0; streams.num_records];
        for record_num in 0..streams.num_records {
            record_num_sorted[record_num] = tuple_attribute_data[record_num].0;
        }

        
        if verbose {
            println!("Reading DEM data...")
        };
        // Load reference raster
        let mut dem = Raster::new(&dem_file, "r")?;

        let rows_dem = dem.configs.rows as isize;
        let columns_dem = dem.configs.columns as isize;

        // Set intermediary arrays and initialize some values
        if verbose {
            println!("Creating intermediary arrays...")
        };
        let mut cost: Array2D<u8> = Array2D::new(rows_dem, columns_dem, 0, 0)?;
        let mut solved: Array2D<u8> = Array2D::new(rows_dem, columns_dem, 0, 0)?;
        let mut accum: Array2D<f64> = Array2D::new(rows_dem, columns_dem, i32::max_value() as f64, -1f64)?;
        let mut backlink: Array2D<u8> = Array2D::new(rows_dem, columns_dem, 255, 255)?;

        let mut bb = BoundingBox {..Default::default()};
        let (mut top_row, mut bottom_row, mut left_col, mut right_col): (isize, isize, isize, isize);
        let mut row_y_coord: f64;
        let mut col_x_coord: f64;
        let (mut x1, mut x2, mut y1, mut y2): (f64, f64, f64, f64);
        let (mut x_prime, mut y_prime): (f64, f64);
        let start_point = 0 as usize;
        let mut end_point: usize;
        let mut output_something = false;
        let num_records = streams.num_records;


        // Set the small increment value to use
        let small_num = if !flat_increment.is_nan() || flat_increment == 0f64 {
            flat_increment
        } else {
            dem.configs.data_type = DataType::F64; // Don't take any chances and promote to 64-bit
            let elev_digits = (dem.configs.maximum as i64).to_string().len();
            let elev_multiplier = 10.0_f64.powi((6 - elev_digits) as i32);
            let resx = dem.configs.resolution_x;
            let resy = dem.configs.resolution_y;
            let diagres = (resx * resx + resy * resy).sqrt();
            1.0_f64 / elev_multiplier as f64 * diagres.ceil()
        };



        // Loop over all stream lines
        // Could be parallelized by chunk of lines of the same priority (or lack of priority
        // if --order_by is not set) but it may not increase speed by much but might worth a try
        for record_num in record_num_sorted {
            let record = streams.get_record(record_num);

            // Ensure that only singlepart polylines are used. This case should be catched right at loading
            // but WhiteboxTools doesn't seem to currently distinguish singlepart from multipart
            // (see lines 500 and up in whitebox-vector\src\shapefile\geometry.rs)
            if record.num_parts as usize > 1
            {
                return Err(Error::new(
                    ErrorKind::InvalidInput,
                    "The input vector data must be of polyline singlepart type, not polyline multipart.",
                ));
            }



            // STREAM RASTERIZATION  (see vector_lines_to_raster.rs for reference)

            end_point = record.num_points as usize - 1;

            bb.initialize_to_inf();
            for i in start_point..end_point + 1 {
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
            if top_row >= rows_dem {
                top_row = rows_dem - 1;
            }
            if bottom_row >= rows_dem {
                bottom_row = rows_dem - 1;
            }

            if left_col < 0 {
                left_col = 0;
            }
            if right_col < 0 {
                right_col = 0;
            }
            if left_col >= columns_dem {
                left_col = columns_dem - 1;
            }
            if right_col >= columns_dem {
                right_col = columns_dem - 1;
            }


            // Find each intersection with a row
            for row in top_row..bottom_row + 1 {
                row_y_coord = dem.get_y_from_row(row);
                // find the x-coordinates of each of the line segments
                // that intersect this row's y coordinate
                for i in start_point..end_point {
                    if is_between(row_y_coord, record.points[i].y, record.points[i + 1].y) {
                        y1 = record.points[i].y;
                        y2 = record.points[i + 1].y;
                        if y2 != y1 {
                            x1 = record.points[i].x;
                            x2 = record.points[i + 1].x;

                            // calculate the intersection point
                            x_prime = x1 + (row_y_coord - y1) / (y2 - y1) * (x2 - x1);
                            let col = dem.get_column_from_x(x_prime);

                            cost.set_value(row, col, 1);
                            output_something = true;
                        }

                    }
                }
            }

            // Find each intersection with a column
            for col in left_col..right_col + 1 {
                col_x_coord = dem.get_x_from_column(col);
                // find the y-coordinates of each of the line segments
                // that intersect this column's x coordinate
                for i in start_point..end_point {
                    if is_between(col_x_coord, record.points[i].x, record.points[i + 1].x) {
                        x1 = record.points[i].x;
                        x2 = record.points[i + 1].x;
                        if x1 != x2 {
                            y1 = record.points[i].y;
                            y2 = record.points[i + 1].y;

                            // calculate the intersection point
                            y_prime = y1 + (col_x_coord - x1) / (x2 - x1) * (y2 - y1);
                            let row = dem.get_row_from_y(y_prime);

                            cost.set_value(row, col, 1);
                            output_something = true;
                        }
                    }
                }
            }

            // Find the row/col intersection point of both ends of the stream
            let row_start = dem.get_row_from_y(record.points[start_point].y);
            let col_start = dem.get_column_from_x(record.points[start_point].x);
            let row_end = dem.get_row_from_y(record.points[end_point].y);
            let col_end = dem.get_column_from_x(record.points[end_point].x);



            // COST DISTANCE ANALYSIS TO FIND LEAST COST PATH  (see cost_distance.rs for reference)

            // Make sure that both ends of the stream are included in the cost raster
            // This is important if two independent streams are connected as the continuity is enforced.
            cost.set_value(row_start, col_start, 1); 
            cost.set_value(row_end, col_end, 1);

            // Set to zero the end of the stream for both the cost accumulation raster
            // and the backlink raster (D8 pointer)
            accum.set_value(row_end, col_end, 0.0);  
            backlink.set_value(row_end, col_end, 0);

            // Cost distance analysis
            let num_cells = ((bottom_row - top_row + 1) * (right_col - left_col + 1)) as usize;
            let mut minheap = BinaryHeap::with_capacity(num_cells);
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
                0f64,
                diag_cell_size,
                cell_size_x,
                diag_cell_size,
                cell_size_y,
                diag_cell_size,
                cell_size_x,
                diag_cell_size,
                cell_size_y,
            ];
            let dx = [0,  1, 1, 1, 0, -1, -1, -1,  0];
            let dy = [0, -1, 0, 1, 1,  1,  0, -1, -1];
            let backlink_dir = [0, 5, 6, 7, 8, 1, 2, 3, 4];
            while !minheap.is_empty() {
                let cell = minheap.pop().expect("Error during pop operation.");
                row = cell.row;
                col = cell.column;
                if solved.get_value(row, col) == 0 {
                    if dem.get_value(row, col) == dem.configs.nodata {
                        return Err(Error::new(
                            ErrorKind::InvalidInput,
                            "The input vector data must not overlap NoData cells from DEM or extend beyond its extent.",
                        ));
                    }
                    solved.set_value(row, col, 1);
                    accum_val = accum.get_value(row, col);
                    cost1 = cost.get_value(row, col).into();
                    for n in 1..9 {
                        col_n = col + dx[n];
                        row_n = row + dy[n];
                        if cost.get_value(row_n, col_n) != cost.nodata {
                            cost2 = cost.get_value(row_n, col_n).into();
                            new_cost = accum_val + (cost1 + cost2) / 2.0 * dist[n];
                            if new_cost < accum.get_value(row_n, col_n) {
                                if solved.get_value(row_n, col_n) == 0 {
                                    accum.set_value(row_n, col_n, new_cost);
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
            


            // DECREASE ELEVATION ALONG BACKLINK (LEAST COST PATH)  (see breach_depressions.rs for reference)

            let mut z_target: f64;
            let mut zin_n: f64;
            let (mut row_n, mut col_n): (isize, isize);
            let mut dir = 1u8;
            row = row_start;
            col = col_start;

            while dir != 0 {
                z_target = dem.get_value(row, col);
                
                // Compare elevation of each neighbour to the central cell
                // and lower it to the bare minimum to clear the lower cells
                for n in 1..9 {
                    row_n = row + dy[n];
                    col_n = col + dx[n];
                    zin_n = dem.get_value(row_n, col_n);
                    if zin_n != dem.configs.nodata {
                        if zin_n <= z_target {
                            z_target = zin_n - small_num;
                        }
                    }
                }
                dem.set_value(row, col, z_target);

                // Get the next cell in line
                dir = backlink.get_value(row, col);
                row += dy[dir as usize];
                col += dx[dir as usize];
            }

            if verbose {
                iter_progess += 1;
                progress = (100.0_f64 * iter_progess as f64 / num_records as f64) as usize;
                if progress != old_progress {
                    println!(
                        "Burning stream {} of {}: {}%",
                        iter_progess,
                        num_records,
                        progress
                    );
                    old_progress = progress;
                }
            }

            // Reset intermediary arrays only in the bbox touched by
            // the current stream to save time
            for row in top_row..bottom_row {
                for col in left_col..right_col {
                    cost.set_value(row, col, cost.nodata);
                    solved.set_value(row, col, solved.nodata);
                    accum.set_value(row, col, i32::max_value() as f64);
                    //backlink.set_value(row, col, backlink.nodata);
                }
            }
        }

        
        // Redirect DEM raster to the output file
        dem.file_name = output_file;
        dem.file_mode = "w".to_string();

        let elapsed_time = get_formatted_elapsed_time(start);
        dem.add_metadata_entry(format!(
            "Created by whitebox_tools\' {} tool",
            self.get_tool_name()
        ));
        dem.add_metadata_entry(format!("Elapsed Time (excluding I/O): {}", elapsed_time));


        let _ = match dem.write() {
            Ok(_) => {
                if verbose {
                    println!("Output file written")
                }
            }
            Err(e) => return Err(e),
        };


        if !output_something && verbose {
            println!("Warning: No streams were burned into the raster.");
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


/// FUNCTION FROM vector_lines_to_raster.rs
fn is_between(val: f64, threshold1: f64, threshold2: f64) -> bool {
    if val == threshold1 || val == threshold2 {
        return true;
    }
    if threshold2 > threshold1 {
        return val > threshold1 && val < threshold2;
    }
    val > threshold2 && val < threshold1
}


/// FUNCTIONS FROM cost_distance.rs
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
