use std::env;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, Write};
use std::path::Path;
use log::{LevelFilter};
use log4rs;

const PARAM_FILE_NAME: &str = "param.txt";
const LOG_FILE_NAME: &str = "pp2_pred_db_log.log";
const DATABASE_FOLDER: &str = "database_folder";

fn main() {
    // Setup logger
    let log_config = log4rs::append::file::FileAppender::builder()
        .encoder(Box::new(log4rs::encode::pattern::PatternEncoder::new("{l} - {m}{n}")))
        .build(LOG_FILE_NAME)
        .unwrap();

    let log_config = log4rs::config::Config::builder()
        .appender(log4rs::config::Appender::builder().build("logfile", Box::new(log_config)))
        .build(log4rs::config::Root::builder().appender("logfile").build(LevelFilter::Info))
        .unwrap();

    log4rs::init_config(log_config).unwrap();

    let input_pdb_given = env::args().nth(1).expect("Please provide a PDB ID as an argument");

    // Load parameter file
    let param_file_path = Path::new(PARAM_FILE_NAME);
    let param_file = File::open(&param_file_path).expect("Error opening parameter file");
    let _param_reader = BufReader::new(param_file);

    // Download or load the input PDB file
    let input_pdb_file_path = format!("{}/{}.pdb", DATABASE_FOLDER, input_pdb_given);
    let _input_pdb_file = download_or_load_pdb(&input_pdb_given, &input_pdb_file_path);

}

fn download_or_load_pdb(pdb_id: &str, pdb_file_path: &str) -> File {
    let path = Path::new(pdb_file_path);
    if path.exists() {
        File::open(&path).expect("Error opening PDB file")
    } else {
        let pdb_url = format!("https://files.rcsb.org/download/{}.pdb", pdb_id);
        let response = reqwest::blocking::get(&pdb_url).expect("Error downloading PDB file");
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .expect("Error creating PDB file");

        write!(file, "{}", response.text().unwrap()).expect("Error writing PDB file");

        file
    }
}


/*
#[macro_use]
extern crate lazy_static;

lazy_static! {
    static ref NEW_LIST_OF_UNIQUE_ALIGNMENT: Mutex<Vec<Vec<String>>> = Mutex::new(Vec::new());
}

const TARGETS_DIR_NAME: &str = "targets";

fn main() {

    // Get the input PDB ID from command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: pp2_pred <PDB code>");
        std::process::exit(1);
    }
    let input_pdb_id = &args[1];

    // Load necessary data
    let (pdb_id, receptor_chain, ppii_chain) = load_data_files();

    // Download PDB file if not already present
    download_pdb_if_not_present(input_pdb_id);

    // Load the input PDB file
    let input_pdb_file_path = format!("database_folder/{}.pdb", input_pdb_id);
    let input_pdb_file = match File::open(&input_pdb_file_path) {
        Ok(file) => file,
        Err(err) => {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        }
    };
    let input_pdb_reader = BufReader::new(input_pdb_file);

    // Parse the input PDB file
    let mut input_pdb_model = Model::new(1);
    for line in input_pdb_reader.lines().flatten() {
        let record = match Record::from_pdb_record(&line) {
            Ok(record) => record,
            Err(_) => continue,
        };
        if let Record::Atom(atom) = record {
            input_pdb_model.add_atom(atom);
        }
    }
    let input_pdb_structure = PDB::new(vec![input_pdb_model], pdb_id.to_owned());

    // Save a new PDB file containing only the specified receptor chain
    let receptor_pdb_file_path = format!("{}_receptor.pdb", input_pdb_id);
    save_pdb(&receptor_pdb_file_path, &input_pdb_structure, &receptor_chain).unwrap();

    // Save a new PDB file containing only the specified PPII chain
    let ppii_pdb_file_path = format!("{}_ppii.pdb", input_pdb_id);
    save_pdb_filtered(&ppii_pdb_file_path, &input_pdb_structure, &ppii_chain).unwrap();

    // Add your own values for donor_dict, predicted_alignments, and files_4_click here
    let donor_dict: HashMap<String, String> = HashMap::new();
    let predicted_alignments: Vec<String> = Vec::new();
    let files_4_click: Vec<String> = Vec::new();
    let input_pdb_given = "1CKA";
    let current_working_dir = "."; // Replace with your own working directory
    let current_fragment = "fragment";

    // Call the neighbor_search function       

    neighbour_search(
        &input_pdb_structure,
        &donor_dict,
        &predicted_alignments,
        &files_4_click,
        input_pdb_given,
        current_working_dir,
        current_fragment,
    );



    if !NEW_LIST_OF_UNIQUE_ALIGNMENT.lock().unwrap().is_empty() {
        let new_list_of_unique_alignment = NEW_LIST_OF_UNIQUE_ALIGNMENT.lock().unwrap();
        let best_alignment = new_list_of_unique_alignment.iter()
            .min_by_key(|x| (
                x[x.len() - 2].parse::<f64>().unwrap(),
                -x[x.len() - 1].parse::<f64>().unwrap()
            ))
            .unwrap();
        let best_alignment_query = [
            best_alignment[0].to_owned(),
            best_alignment[4][0..1].to_owned(),
            best_alignment[4][1..2].to_owned(),
            best_alignment[4][2..].to_owned(),
            best_alignment[5].to_owned(),
            best_alignment[6].to_owned(),
            best_alignment[8][best_alignment[8].len() - 4..].to_owned(),
        ];
        println!("{:?}", *new_list_of_unique_alignment);
        if best_alignment[best_alignment.len() - 2].parse::<f64>().unwrap() >= 2.5 {
            println!("No binding site found in the given query structure {}", input_pdb_given);
            std::process::exit(1);
        }

        println!(
            "For query structure {}, predicted binding site details are - Model = {}, Chain = {}, TRP = {}, NBR = {}",
            best_alignment[0],
            best_alignment[4][0..1],
            best_alignment[4][1..2],
            best_alignment[4][2..],
            best_alignment[5]
        );
        println!(
            "Template PPII is {} with a Score of {}",
            best_alignment[8][best_alignment[8].len() - 4..],
            best_alignment[best_alignment.len() - 2]
        );
    } else {
        println!("No binding site found in the given query structure {}", input_pdb_given);
        std::process::exit(1);
    }
}



fn glob(folder: &str, pattern: &str) -> Vec<String> {
    let pattern = format!("{}/{}", folder, pattern);
    let mut result = Vec::new();

    for entry in glob::glob(&pattern)? {
        match entry {
            Ok(path) => result.push(path.to_str().unwrap().to_owned()),
            Err(e) => println!("{:?}", e),
        }
    }
    result
}

fn first<T>(items: &mut Vec<T>) -> Option<T> {
    if items.is_empty() {
        None
    } else {
        Some(items.remove(0))
    }
}


fn click4all(input_pdb1: &str, input_pdb2: &str) {
    let cmd = format!("./click {} {} >/dev/null 2>&1", input_pdb1, input_pdb2);
    let output = Command::new("sh")
        .arg("-c")
        .arg(&cmd)
        .output()
        .expect("failed to execute process");
    assert!(output.status.success());
}


// This function takes the input PDB ID and an input chain, and saves a new PDB file
// containing only the specified chain
fn save_pdb(filename: &str, pdb_id: &str, chain: &str) -> std::io::Result<()> {
    let mut output_file = File::create(filename)?;
    for model in pdb.models() {
        for chain_data in model.chain_iter(chain) {
            let chain_model = Model::new(1);
            for residue in chain_data.residues() {
                for atom in residue.atoms() {
                    chain_model.add_atom(atom.clone());
                }
            }
            chain_model.write(&mut output_file)?;
        }
    }
    Ok(())
}

fn load_data_files() -> (String, String, String) {
    // Load necessary data files
    // Return a tuple with pdb_id, receptor_chain, and ppii_chain
    
    let pdb_id = "1CKA".to_string();
    let receptor_chain = "A".to_string();
    let ppii_chain = "B".to_string();

    (pdb_id, receptor_chain, ppii_chain)
}


fn download_pdb_if_not_present(input_pdb_given: &str) {
    let pdb_file_path = format!("database_folder/{}.pdb", input_pdb_given);
    if !Path::new(&pdb_file_path).exists() {
        let pdb_url = format!("https://files.rcsb.org/download/{}.pdb", input_pdb_given);
        let mut response = reqwest::get(&pdb_url).expect("Failed to download PDB file");
        let mut pdb_file = File::create(&pdb_file_path).expect("Failed to create PDB file");
        io::copy(&mut response, &mut pdb_file).expect("Failed to save PDB file");
    }
}

fn save_pdb(filename: &str, pdb_id: &str, chain: char) -> std::io::Result<()> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let output_file_path = format!("{}_4sim.pdb", pdb_id);
    let mut output_file = File::create(&output_file_path)?;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("ATOM") {
            let chain_name = line.chars().nth(21).unwrap();
            if chain_name == chain {
                writeln!(output_file, "{}", line)?;
            }
        }
    }

    Ok(())
}

const NEIGHBORHOOD_LOOK_UP_CUT_OFF: f64 = 12.0;
const H_BOND_CUT_OFF: f64 = 3.5;

fn calc_angle(a: Array2<f64>, b: Array2<f64>, c: Array2<f64>) -> f64 {
    let ba = &a - &b;
    let bc = &c - &b;
    let cosine_angle = ba.dot(&bc) / (ba.norm() * bc.norm());
    cosine_angle.acos()
}



fn neighbour_search(
    structure: &pdb::PDB,
    donor_dict: &HashMap<(String, String), String>,
    predicted_alignments: &Vec<String>,
    files_4_click: &Vec<String>,
    input_pdb_given: &str,
    current_working_dir: &str,
    current_fragment: &str,
) {
    let model = structure.model(1).unwrap();
    for chain in model.chains() {
        for residue in chain.residues() {
            if residue.resname() == "TRP" {
                let the_ne1_atom = residue.atom("NE1").unwrap();
                println!(
                    "Residue Tryptophan is present at: {}, {}, {}",
                    the_ne1_atom.model(),
                    the_ne1_atom.chain(),
                    the_ne1_atom.residue()
                );

                let chain_atoms: Vec<_> = model.atoms().collect();
                let neighbourhood_search = NeighborSearch::new(&chain_atoms, NEIGHBORHOOD_LOOK_UP_CUT_OFF);
                let neighbour_atoms = neighbourhood_search.neighbors(&the_ne1_atom);

                for n_atom in neighbour_atoms {
                    let mut rejection_list = vec![];

                    let atom_dic = (n_atom.residue().resname().to_string(), n_atom.name().to_string());
                    if n_atom != the_ne1_atom && donor_dict.contains_key(&atom_dic) {
                        let internal_look_up = neighbourhood_search.neighbors_within(&n_atom, H_BOND_CUT_OFF);

                        for internal_atoms in internal_look_up {
                            if internal_atoms != n_atom {
                                let n_atom_residue = n_atom.residue().resname();
                                let internal_look_up_residue = internal_atoms.residue().resname();

                                let n_atom_id = format!("{}:{}", n_atom.name(), n_atom_residue);
                                let internal_look_up_residue_id = format!("{}:{}", internal_atoms.name(), internal_look_up_residue);

                                // Handle the internal hydrogen bond look up here
                                // ...
                            }
                        }

                        if !rejection_list.contains(&n_atom) && n_atom != the_ne1_atom && donor_dict.contains_key(&atom_dic) {
                            // Carve the segment of query pdb id for CLICK alignment.
                            // ...

                            for dataset_file in all_data_files.iter() {
                                if dataset_file.len() == 16 {
                                    let the_path = format!(
                                        "{}_{}_{}_{}_{}_{}",
                                        n_atom.model(),
                                        n_atom.chain(),
                                        n_atom.residue(),
                                        the_ne1_atom.serial(),
                                        n_atom.serial(),
                                        n_atom.residue().resname()
                                    );

                                    // Mask Template PDB atoms for CLICK alignment
                                    // ...

                                    // Mask Query Segment atoms which was "Carved" for CLICK alignment
                                    // ...
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

*/


