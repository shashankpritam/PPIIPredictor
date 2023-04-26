use std::env;
use std::error::Error;
use std::path::Path;
use std::io::{BufReader, BufRead, Read, Write, BufWriter};
use std::fs::{self, File, create_dir_all};
use std::process::Command;
use std::collections::HashMap;
use glob::glob;
use reqwest;
use pp2predictor::pdb_parser::{parse_pdb_file, Structure, Model, Chain, Residue, Atom, write_pdb, NeighborSearch, calculate_distance, calculate_angle};
use std::f32::consts::PI;
use std::rc::Rc;
use std::collections::HashSet;
use lazy_static::lazy_static;

// In src/main.rs
use kdtree::KdTree;
use kdtree::distance::squared_euclidean;
use std::sync::Arc;

lazy_static! {
    static ref ALL_DATASET_PDB: Vec<String> = {
        glob("dataset/*.pdb")
            .expect("Failed to read glob pattern")
            .map(|entry| entry.expect("Failed to read entry").display().to_string())
            .collect()
    };
}

lazy_static! {
    static ref DONOR_DICT: HashSet<(&'static str, &'static str)> = {
        let donors = [
            ("ARG", "NE"),
            ("ARG", "NH1"),
            ("ARG", "NH2"),
            ("ASN", "ND2"),
            ("ASX", "ND2"),
            ("CYS", "SG"),
            ("GLN", "NE2"),
            ("GLX", "NE2"),
            ("HIS", "ND1"),
            ("HSE", "NE2"),
            ("HSP", "ND1"),
            ("HSP", "NE2"),
            ("LYS", "NZ"),
            ("SER", "OG"),
            ("THR", "OG1"),
            ("TRP", "NE1"),
            ("TYR", "OH"),
        ];
        donors.iter().cloned().collect()
    };
}


// Add undeclared constant values
const LOG_FILE_NAME: &str = "log.txt";
const PARAM_FILE_NAME: &str = "param.txt";
const DATABASE_FOLDER: &str = "pdb_files";
const H_BOND_ACCEPTOR_LIST: &str = "H_bond_Acceptor_List";
const H_BOND_DONOR_LIST: &str = "H_bond_Donor_List";
const H_BOND_INFO: &str = "data_hbond/hbond_trp_all.txt";

const PDB_INFO: &[(&str, &str, &str)] = &[
    ("4B6H", "A", "C"), ("4WSF", "A", "B"), ("5J3T", "A", "C"), 
    ("1L2Z", "A", "B"), ("3FMA", "A", "L"), ("4CW1", "A", "F"), 
    ("5C0D", "A", "C"), ("1JK8", "B", "C"), ("3PDO", "B", "C"), 
    ("4P57", "A", "B"), ("2ODD", "A", "B"), ("2V8F", "B", "C"), 
    ("1CKA", "A", "B"), ("1GBQ", "A", "B"), ("1SEM", "A", "C"), 
    ("1UTI", "A", "D"), ("1YWO", "A", "P"), ("2DF6", "A", "C"), 
    ("2DRM", "A", "E"), ("2J6F", "A", "C"), ("2LCS", "A", "B"), 
    ("2ROL", "A", "B"), ("2RPN", "A", "B"), ("2VKN", "A", "C"), 
    ("2VWF", "A", "B"), ("3I5R", "A", "B"), ("3U23", "A", "B"), 
    ("3ULR", "B", "C"), ("4CC2", "A", "B"), ("4F14", "A", "B"), 
    ("4HVW", "A", "B"), ("4J9C", "A", "B"), ("4LNP", "A", "B"), 
    ("4U5W", "D", "C"), ("1JMQ", "A", "P"), ("2EZ5", "W", "P"), 
    ("2JO9", "A", "B"), ("2LAW", "A", "B"), ("2LAZ", "A", "B")
];



fn main() {
    let args: Vec<String> = env::args().collect();
    let input_pdb_given = args
        .get(1)
        .expect("Please provide a PDB ID as an argument")
        .to_owned();

    // Add model id as input with default value 0
    let model_id: usize = args
        .get(2)
        .map_or_else(|| Ok(0), |arg| usize::from_str_radix(arg, 10))
        .unwrap_or(0);

    let param_file_path = Path::new(PARAM_FILE_NAME);
    let param_file = File::open(&param_file_path).expect("Error opening parameter file");
    let _param_reader = BufReader::new(param_file);


    let input_pdb_file_path = format!("{}/{}.pdb", DATABASE_FOLDER, input_pdb_given);
    let _input_pdb_file = download_or_load_pdb(&input_pdb_given, &input_pdb_file_path);

    // Load PDB chains
    let (pdb_id, _receptor_chain, _ppii_chain) = load_pdb_chains();

    // Get peptide chain and receptor chain of input PDB
    let input_pdb_upper = input_pdb_given.to_uppercase();
    let pdb_id_index = pdb_id.iter().position(|id| id == &input_pdb_upper).unwrap();
    let _input_pdb = &pdb_id[pdb_id_index];
    let _pep_chain = get_pep_chain(&input_pdb_given);
    let _rec_chain = get_rec_chain(&input_pdb_given);

    // Read the content of the PDB file
    let file_content = fs::read_to_string(&input_pdb_file_path).expect("Error reading PDB file");

    // Parse the content of the PDB file
    let query_structures =
        parse_pdb_file(&file_content).expect("Error parsing PDB file");
    //println!("Structures: {:?}", query_structures);
    let neighbourhood_look_up_cut_off = 3.5;
    let h_bond_cut_off = 2.5;
    let all_dataset_pdb: Vec<&str> = ALL_DATASET_PDB.iter().map(|s| s.as_ref()).collect();
    neighbour_search(&query_structures, &all_dataset_pdb);
}



fn download_or_load_pdb(pdb_id: &str, pdb_file_path: &str) -> File {
    let path = Path::new(pdb_file_path);
    let alternative_pdb_file_path = format!("{}/{}.ent", DATABASE_FOLDER, pdb_id);
    let alternative_path = Path::new(&alternative_pdb_file_path);

    if path.exists() {
        println!("Found PDB file at path: {:?}", &path);
        File::open(&path).expect("Error opening PDB file")
    } else if alternative_path.exists() {
        println!("Found alternative PDB file at path: {:?}", &alternative_path);
        File::open(&alternative_path).expect("Error opening alternative PDB file")
    } else {
        println!("Downloading PDB file from URL: https://files.rcsb.org/download/{}.pdb", pdb_id);
        let pdb_url = format!("https://files.rcsb.org/download/{}.pdb", pdb_id);
        let response = reqwest::blocking::get(&pdb_url).expect("Error downloading PDB file");
        let mut file = std::fs::File::create(&path).expect("Error creating PDB file");

        file.write_all(response.text().unwrap().as_bytes())
            .expect("Error writing PDB file");
        println!("PDB file downloaded and saved at path: {:?}", &path);
        file
    }
}


// A function to get the peptide chain of a PDB from the template dataset containing 39 PDBs
fn get_pep_chain(input_pdb: &str) -> String {
    let pdb_id_upper = input_pdb.to_uppercase();
    let input_peptide_chain = PDB_INFO
        .iter()
        .find(|&&(id, _, _)| id == pdb_id_upper)
        .map(|&(_, _, pep_chain)| pep_chain.to_string())
        .unwrap();
    input_peptide_chain
}

// A function to get the receptor chain of a PDB from the template dataset containing 39 PDBs
fn get_rec_chain(input_pdb: &str) -> String {
    let pdb_id_upper = input_pdb.to_uppercase();
    let input_receptor_chain = PDB_INFO
        .iter()
        .find(|&&(id, _, _)| id == pdb_id_upper)
        .map(|&(_, rec_chain, _)| rec_chain.to_string())
        .unwrap();
    input_receptor_chain
}

fn load_pdb_chains() -> (Vec<String>, Vec<String>, Vec<String>) {
    let (pdb_id, receptor_chain, ppii_chain) = PDB_INFO.iter().fold(
        (Vec::new(), Vec::new(), Vec::new()),
        |(mut pdb_id, mut receptor_chain, mut ppii_chain), &(id, rec_chain, pep_chain)| {
            pdb_id.push(id.to_string());
            receptor_chain.push(rec_chain.to_string());
            ppii_chain.push(pep_chain.to_string());
            (pdb_id, receptor_chain, ppii_chain)
        },
    );

    (pdb_id, receptor_chain, ppii_chain)
}



fn read_acceptor_list(content: &str) -> HashMap<String, (String, String)> {
    let mut acceptor_list = HashMap::new();

    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        acceptor_list.insert(parts[0].to_string(), (parts[1].to_string(), parts[2].to_string()));
    }

    acceptor_list
}

fn read_donor_list(content: &str) -> HashMap<String, String> {
    let mut donor_list = HashMap::new();

    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        donor_list.insert(parts[0].to_string(), parts[1].to_string());
    }

    donor_list
}




// This function runs the CLICK program on two input PDB files and discards the output
fn click4all(input_pdb1: &str, input_pdb2: &str) {
    let cmd = format!("./click {} {} > /dev/null 2>&1", input_pdb1, input_pdb2);
    Command::new("bash")
        .arg("-c")
        .arg(&cmd)
        .output()
        .expect("Error running CLICK command.");
}

// This function takes an input PDB ID from the template dataset and returns
// the list of TRP involved in Hydrogen Bond with the PPII - Only for the 39 PDB in dataset
// See comment above hbond_files
fn hbond_trp(input_pdb: &str) -> Vec<String> {
    let mut list_of_trp: Vec<String> = Vec::new();
    if let Ok(file) = File::open(H_BOND_INFO) {
        let lines = BufReader::new(file).lines().filter_map(|line| line.ok());

        for line in lines {
            let data: Vec<&str> = line.split_whitespace().collect();
            if input_pdb == data[0] && data.contains(&"TRP") {
                if let Some(trp_res_id) = data
                    .iter()
                    .position(|&s| s == "TRP")
                    .and_then(|pos| data.get(pos + 1))
                {
                    let trp_res_id_string = trp_res_id.to_string();
                    if !list_of_trp.contains(&trp_res_id_string) {
                        list_of_trp.push(trp_res_id_string);
                    }
                }
            }
        }
    } else {
        panic!("Error opening hbond_trp_all.txt file");
    }
    list_of_trp
}


// This function takes input of template pdb id and an side chain donor atom (scda), side chain donor residue (scdr) of the scda
// and the suffix which acts as identifier as template file. The last parameter save_path is path where the ouput file is
// supposed to be saved.
// This function replaces the scda of scdr with NX and CB atom of scdr as EE.
// This function also replaces the NE1, CA and CZ3 of the Tryptophan residues (TRPs) from the trp_list acquired from hbond_trp
// as AA, BB and CC respectively.
// After replacing the atom names this functions saves the new PDB as per provided inputs.
// Suffix is used for renaming the output file.
// SCDA = Side Chain Donor ATOM
// SCDR = Side Chain Donor RESIDUE
fn mask_temp_atoms(structure: &mut Structure, scda: &str, scdr: &str, suffix: &str, save_path: &str, the_pdb_id: &str,) -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    let input_pdb_given = args
        .get(1)
        .expect("Please provide a PDB ID as an argument")
        .to_owned();
    let save_dir = Path::new("click_output").join(save_path);
    create_dir_all(&save_dir)?;
    let file_name = format!("{}{}{}{}{}{}.pdb", input_pdb_given, the_pdb_id, scda, scdr, save_path, suffix);
    let output_pdb_path = save_dir.join(file_name);

    let mut output_pdb_file = File::create(output_pdb_path)?;
    for model in &mut structure.models {
        for chain in &mut model.chains {
            for residue in &mut chain.residues {
                let res_name = &residue.name;
                let res_seq = residue.id;
                for atom in &mut residue.atoms {
                    let atm_name = &atom.name;
                    match (atm_name.as_str(), res_name.as_str()) {
                        ("NE1", "TRP") => atom.name = "AA ".to_string(),
                        ("CA", "TRP") => atom.name = "BB".to_string(),
                        ("CZ3", "TRP") => atom.name = "CC ".to_string(),
                        _ if atm_name == scda && res_name == scdr => match atm_name.len() {
                            1 => atom.name = "NX".to_string(),
                            2 => atom.name = "NX".to_string(),
                            3 => atom.name = "NX ".to_string(),
                            _ => atom.name = " NX ".to_string(),
                        },
                        _ if atm_name == "CB" && res_name == scdr => atom.name = "EE".to_string(),
                        _ => (),
                    }
                }
            }
        }
    }

    // Write the modified structure to the output file
    write!(output_pdb_file, "{:#?}", structure)?;

    Ok(())
}


// This function takes input of query pdb id and an Tryptophan residue as Biopython Residue Object - the_trp
// Other input parameters are - the_nbr; which the scdr but as biopython RESIDUE object and the_nbr_dnr -
// which the scda but as biopython ATOM object
// The suffix which acts as identifier as query file. The last parameter save_path is path where the ouput file is
// supposed to be saved.
// This function replaces the the_nbr with NX and CB atom of the_nbr_dnr as EE.
// This function also replaces the NE1, CA and CZ3 of the Tryptophan residues (TRPs) from the_trp
// as AA, BB and CC respectively.
// After replacing the atom names this functions saves the new PDB as per provided inputs.
// Suffix is used for renaming the output file.
fn mask_query_atoms(input_pdb: &str, the_trp: &str, the_nbr: &str, the_nbr_dnr: &str, suffix: &str, save_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let save_dir = Path::new("click_output").join(save_path);
    create_dir_all(&save_dir)?;

    let file_name = format!("{}{}.pdb", &input_pdb[..input_pdb.len() - 4], suffix);
    let output_pdb_path = save_dir.join(file_name);

    let input_pdb_file = File::open(input_pdb)?;
    let reader = BufReader::new(input_pdb_file);

    let mut output_pdb_file = File::create(output_pdb_path)?;

    let query_structure = parse_pdb_file(input_pdb)?;

    let mut the_trp_residue = None;
    let mut the_nbr_residue = None;

    for model in &query_structure.models {
        for chain in &model.chains {
            for residue in &chain.residues {
                if residue.name == the_trp {
                    the_trp_residue = Some(residue);
                } else if residue.name == the_nbr {
                    the_nbr_residue = Some(residue);
                }
            }
        }
    }

    let the_trp_residue = the_trp_residue.unwrap();
    let the_nbr_residue = the_nbr_residue.unwrap();
    let the_nbr_residue_int = the_nbr_residue.id as i32;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("ATOM") {
            let atm_name = line[12..16].trim();
            let res_seq = line[22..26].trim();

            let mut new_line = line.clone();
            if atm_name == "NE1" && res_seq.parse::<i32>()? == the_trp_residue.id as i32 {
                new_line = line.replace(atm_name, "AA ");
            } else if atm_name == "CA" && res_seq.parse::<i32>()? == the_trp_residue.id as i32 {
                new_line = line.replace(atm_name, "BB");
            } else if atm_name == "CZ3" && res_seq.parse::<i32>()? == the_trp_residue.id as i32 {
                new_line = line.replace(atm_name, "CC ");
            } else if atm_name == the_nbr_dnr && res_seq.parse::<i32>()? == the_nbr_residue_int {
                match atm_name.len() {
                    1 => new_line = line.replacen(atm_name, "NX", 1),
                    2 => new_line = line.replacen(atm_name, "NX", 1),
                    3 => new_line = line.replacen(atm_name, "NX ", 1),
                    _ => new_line = line.replacen(atm_name, " NX ", 1),
                }
            } else if atm_name == "CB" && res_seq.parse::<i32>()? == the_nbr_residue_int {
                new_line = line.replace(atm_name, "EE");
            }

            writeln!(output_pdb_file, "{}", new_line)?;
        } else {
            writeln!(output_pdb_file, "{}", line)?;
                    }
    }

    Ok(())
}

// This function "undoes" what masking function above does.
// That is after the CLICK alignments all the desired files will have "normal"/unmasked file n_atom_residue
// with _new as suffix for identifier.
fn unmask_atoms_save(input_pdb: &str, input_chain: char) -> Result<(), Box<dyn std::error::Error>> {
    let input_pdb_file = File::open(input_pdb)?;
    let reader = BufReader::new(input_pdb_file);

    let output_pdb_path = format!("{}_new.pdb", &input_pdb[..input_pdb.len() - 4]);
    let mut output_pdb_file = File::create(&output_pdb_path)?;

    let file_path = Path::new(input_pdb).canonicalize()?;
    let the_dn = file_path.file_stem().unwrap().to_str().unwrap().split('_').last().unwrap();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("ATOM") {
            let atm_name = line[12..16].trim();
            let chain_id = line[21..22].trim().chars().next().unwrap();

            let mut new_line = line.clone();
            if chain_id == input_chain {
                match atm_name {
                    "AA" => new_line = line.replace("AA", "NE1"),
                    "BB" => new_line = line.replace("BB", "CA"),
                    "CC" => new_line = line.replace("CC", "CZ3"),
                    "NX" => new_line = line.replacen("NX", the_dn, 1),
                    "EE" => new_line = line.replace("EE", "CB"),
                    _ => (),
                }
            }

            writeln!(output_pdb_file, "{}", new_line)?;
        } else {
            writeln!(output_pdb_file, "{}", line)?;
        }
    }

    Ok(())
}

// This function takes a pdb filename with complete path, it is 4 letter ID and a Chain ID as as input.
// This function then returns a new file in the parent directory with _4sim as suffix identifier.
// That PDB file should contain the pdb id provided with only the chain which was provided as an input parameter.
fn save_pdb(filename: &str, pdb_id: &str, chain: char) -> Result<(), Box<dyn std::error::Error>> {
    let input_pdb_file = File::open(filename)?;
    let reader = BufReader::new(input_pdb_file);

    let current_working_dir = std::env::current_dir()?;
    let output_pdb_path = current_working_dir.join(format!("{}_4sim.pdb", pdb_id));
    let mut output_pdb_file = File::create(output_pdb_path)?;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with("ATOM") {
            let chain_name = line[21..22].trim().chars().next().unwrap();

            if chain_name == chain {
                writeln!(output_pdb_file, "{}", line)?;
            }
        }
    }

    Ok(())
}



pub fn carve(
    input_file_path: &str,
    output_file_path: &str,
    chain_ids: &[char],
    residue_ids: &[i32],
) -> Result<(), Box<dyn Error>> {
    // Open input file
    let input_file = File::open(input_file_path)?;
    let mut input_content = String::new();
    let mut reader = BufReader::new(input_file);
    reader.read_to_string(&mut input_content)?;
    let structure = parse_pdb_file(&input_content)?;

    // Filter the structure based on the given chain IDs and residue IDs
    let mut filtered_structure = Structure {
        models: Vec::new(),
        chains: Vec::new(),
    };
    for model in structure.models {
        let mut chains = Vec::new();
        for chain in model.chains {
            if chain_ids.contains(&chain.id) {
                let mut residues = Vec::new();
                for residue in chain.residues {
                    if residue_ids.contains(&(residue.id as i32)) {
                        residues.push(residue);
                    }
                }
                if !residues.is_empty() {
                    chains.push(Chain {
                        id: chain.id,
                        residues,
                    });
                }
            }
        }
        if !chains.is_empty() {
            filtered_structure.models.push(Model {
                chains,
                ..model
            });
        }
    }

    // Open output file
    let output_file = File::create(output_file_path)?;
    let mut output_pdb_file = BufWriter::new(output_file);

    // Write filtered structure to output file
    write_pdb(&filtered_structure, &mut output_pdb_file)?;

    Ok(())
}





// Define a custom atom structure to hold a reference to its parent residue
struct AtomWithParent<'a> {
    atom: &'a Atom,
    parent_residue: &'a Residue,
}


fn neighbour_search(structure: &Structure, all_dataset_pdb: &[&str]) -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    let input_pdb_given = args
        .get(1)
        .expect("Please provide a PDB ID as an argument")
        .to_owned();

    let mut all_atoms: Vec<Atom> = Vec::new();
    for model in &structure.models {
        for chain in &model.chains {
            for residue in &chain.residues {
                all_atoms.extend(residue.atoms.clone());
            }
        }
    }

    let neighbor_search = NeighborSearch::new(&all_atoms);

    for model in structure.models.iter() {
        for chain in model.chains.iter() {
            for residue in chain.residues.iter() {
                if residue.name == "TRP" {
                    println!("Found Tryptophan residue");
                    if let Some(the_ne1_atom) = residue.get_atom_by_name("NE1") {
                        println!(
                            "Residue Tryptophan is present at: chain {}, residue name {}, atom serial {}",
                            the_ne1_atom.chain_id, the_ne1_atom.res_name, the_ne1_atom.serial
                        );

                        let neighbors = neighbor_search.search_neighbors(&the_ne1_atom, 12.0);
                        let mut unique_neighbors: HashSet<(String, String, char, isize, isize, isize)> = HashSet::new();
                        let mut unique_donor_neighbors: HashSet<(String, String, char, isize, isize, isize)> = HashSet::new();

                        for neighbor in neighbors {
                            let distance = calculate_distance(&the_ne1_atom, &neighbor);

                            if distance > 0.0 {
                                let neighbor_tuple = (
                                    neighbor.name.clone(),
                                    neighbor.res_name.clone(),
                                    neighbor.chain_id,
                                    model.serial_number, // Corrected model serial number field
                                    neighbor.res_seq, // Corrected residue ID field
                                    neighbor.serial, // Added atom serial number
                                );

                                if DONOR_DICT.contains(&(neighbor.res_name.as_str(), neighbor.name.as_str())) {
                                    unique_donor_neighbors.insert(neighbor_tuple.clone());
                                }

                                unique_neighbors.insert(neighbor_tuple);
                            }
                        }

                        for (name, res_name, chain_id, model_serial, residue_id, atom_serial) in unique_donor_neighbors {
                            println!(
                                "Neighbor found: Model {}, Atom {}, Residue {}, Chain {}, Residue ID {}, Atom Serial {}",
                                model_serial, name, res_name, chain_id, residue_id, atom_serial
                            );
                            for dataset_file in ALL_DATASET_PDB.iter() {
                                let the_pdb_id = &dataset_file[(dataset_file.len() - 8)..(dataset_file.len() - 4)];
                                let file_name = format!(
                                    "{}_{}_{}_{}_{}_{}_{}_{}_{}_{}",
                                    input_pdb_given,
                                    model.serial_number,
                                    the_ne1_atom.chain_id,
                                    residue.id,
                                    the_ne1_atom.res_name,
                                    "NE1",
                                    chain_id,
                                    residue_id,
                                    res_name,
                                    name,
                                );
                                println!("{}", file_name);


                                let file_content_to_modify = fs::read_to_string(&dataset_file).expect("Error reading PDB file");
                                let mut structure_to_modify =
                                    parse_pdb_file(&file_content_to_modify).expect("Error parsing PDB file");
                                //mask_temp_atoms(&mut structure_to_modify, &name, &res_name, "_rnmd_ds", &the_path, &the_pdb_id)?;
                            }
                        }


                    } else {
                        println!("NE1 atom not found in Tryptophan residue");
                    }
                }
            }
        }
    }
    Ok(())
}


