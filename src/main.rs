use std::env;
use std::error::Error;
use std::path::Path;
use std::io::{BufRead, BufReader};
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use reqwest;
use std::convert::TryInto;
use pp2predictor::pdb_parser::{parse_pdb_file, Structure};


const DATABASE_FOLDER: &str = "pdb_files";
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
    let input_pdb_given = env::args()
        .nth(1)
        .expect("Please provide a PDB ID as an argument");

    let input_pdb_file_path = format!("{}/{}.pdb", DATABASE_FOLDER, input_pdb_given);
    let _input_pdb_file = download_or_load_pdb(&input_pdb_given, &input_pdb_file_path);

    // Read the content of the PDB file
    let file_content = fs::read_to_string(&input_pdb_file_path).expect("Error reading PDB file");

    // Parse the content of the PDB file
    let query_structures =
        parse_pdb_file(&file_content).expect("Error parsing PDB file");
    println!("Structures: {:?}", query_structures);

    for (model_index, model) in query_structures.models.iter().enumerate() {
        //println!("Model {}: {:?}", model_index + 1, model);
        for (chain_index, chain) in model.chains.iter().enumerate() {
            //println!("  Chain {}: {:?}", chain_index + 1, chain);
            for (residue_index, residue) in chain.residues.iter().enumerate() {
                //println!("    Residue {}: {:?}", residue_index + 1, residue);
                for (atom_index, atom) in residue.atoms.iter().enumerate() {
                    println!("      Atom {}: {:?}", atom_index + 1, atom);
                     println!(" Booyah");
                }
            }
        }
    }
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
