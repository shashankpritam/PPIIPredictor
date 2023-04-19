use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;
use pdb::{Atom, Model, PDB, Record};
use std::process::Command;
use std::path::Path;



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
    for line in input_pdb_reader.lines() {
        let line = match line {
            Ok(line) => line,
            Err(err) => {
                eprintln!("Error: {}", err);
                std::process::exit(1);
            }
        };
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
    save_pdb(&ppii_pdb_file_path, &input_pdb_structure, &ppii_chain).unwrap();
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
fn save_pdb(filename: &str, pdb: &PDB, chain: &str) -> std::io::Result<()> {
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

fn neighbor_search(
    structure: &PDB,
    donor_dict: &HashMap<String, String>,
    predicted_alignments: &Vec<String>,
    files_4_click: &Vec<String>,
    input_pdb_given: &str,
    current_working_dir: &str,
    current_fragment: &str,
) {
    for (model_idx, model) in structure.models.iter().enumerate().take(1) {
        for (chain_idx, chain) in model.iter().enumerate() {
            for (residue_idx, residue) in chain.iter().enumerate() {
                if residue.name() != "TRP" {
                    continue;
                }

                let the_ne1_atom = match residue.get("NE1") {
                    Some(atom) => atom,
                    None => continue,
                };
                println!(
                    "Residue Tryptophan is present at : {:?} {} {} {}",
                    structure.path, model_idx, chain_idx, residue_idx
                );

                let chain_atoms = model.atoms();
                let mut neighborhood_search = NeighborSearch::new(chain_atoms);
                let neighbor_atoms = neighborhood_search.search(the_ne1_atom, NEIGHBORHOOD_LOOK_UP_CUT_OFF);

                for n_atom in neighbor_atoms {
                    let atom_dic = (n_atom.residue().name().to_string(), n_atom.name().to_string());
                    if n_atom == the_ne1_atom || !donor_dict.contains_key(&atom_dic.0) {
                        continue;
                    }

                    let internal_lookup = neighborhood_search.search(n_atom, H_BOND_CUT_OFF);
                    for internal_atoms in internal_lookup {
                        if internal_atoms == n_atom {
                            continue;
                        }

                        let n_atom_residue = n_atom.residue().name();
                        let internal_lookup_residue = internal_atoms.residue().name();
                        let n_atom_id = format!("{}:{}", n_atom.name(), n_atom_residue);
                        let internal_lookup_residue_id = format!("{}:{}", internal_atoms.name(), internal_lookup_residue);

                        let list_of_unique_alignment = predicted_alignments.iter()
                            .filter(|alignment| {
                                let alignments: Vec<&str> = alignment.split('_').collect();
                                alignments[0].to_lowercase() != alignments[8][..4].to_lowercase()
                            })
                            .collect::<Vec<_>>();

                        let current_input_pdb = input_pdb_given.to_lowercase();
                        let click_output_dir = format!("{}/click_output/", current_working_dir);
                        for folder in files_4_click {
                            if !folder.to_lowercase().starts_with(&click_output_dir) {
                                continue;
                            }

                            let dataset_renamed_file = glob(folder, "*_rnmd.1.pdb");
                            let renamed_pdb = glob(folder, "*_rnmd_ds.1.pdb");
                            let click_file = glob(folder, "*.clique");

                            let carved_frag_info: Vec<String> = click_file[0]
                                .split('/')
                                .last()
                                .unwrap()
                                .split("_")
                                .map(|s| s.to_owned())
                                .collect();

                            let carved_frag_info = vec![                                carved_frag_info[0].to_owned(),
                                carved_frag_info[4][0..1].to_owned(),
                                carved_frag_info[4][1..2].to_owned(),
                                carved_frag_info[4][2..].to_owned(),
                                carved_frag_info[5].to_owned(),
                                carved_frag_info[6].to_owned(),
                                carved_frag_info[8][5..9].to_owned(),
                            ];

                            let dataset_rnmd_file = dataset_renamed_file.first();
                            let rnmd_ds = renamed_pdb.first()
                            let dataset_rnmd_file = dataset_renamed_file.first();
                            let renamed_pdb_file = renamed_pdb.first();
                            let click_file = click_files.first();

                            if let (Some(dataset_rnmd_file), Some(renamed_pdb_file), Some(click_file)) = (dataset_rnmd_file, renamed_pdb_file, click_file) {
                                let carved_frag_info: Vec<String> = click_file
                                    .split('/')
                                    .last()
                                    .unwrap()
                                    .split("_")
                                    .map(|s| s.to_owned())
                                    .collect();

                                let carved_frag_info = vec![        carved_frag_info[0].to_owned(),
                                    carved_frag_info[4][0..1].to_owned(),
                                    carved_frag_info[4][1..2].to_owned(),
                                    carved_frag_info[4][2..].to_owned(),
                                    carved_frag_info[5].to_owned(),
                                    carved_frag_info[6].to_owned(),
                                    carved_frag_info[8][5..9].to_owned(),
                                ];

                                let rnmd_ds = renamed_pdb_file.to_str().unwrap();
                                let mut node_id = carved_frag_info[6].to_owned();
                                let full_id = format!(
                                    "{}_{}_{}_{}",
                                    carved_frag_info[1], carved_frag_info[2], carved_frag_info[3], carved_frag_info[5]
                                );
                                let target_file_name = format!(
                                    "{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}",
                                    rnmd_ds.split("/").last().unwrap().split("_").nth(0).unwrap(),
                                    full_id,
                                    input_pdb_given,
                                    current_fragment,
                                    carved_frag_info[0],
                                    carved_frag_info[1],
                                    carved_frag_info[2],
                                    carved_frag_info[3],
                                    carved_frag_info[5],
                                    node_id,
                                    align_query[0],
                                    align_query[6],
                                    align_query[4]
                                );

                                let target_file_path = format!("{}/{}/", current_working_dir.clone(), TARGETS_DIR_NAME);
                                let mut full_target_file_path = target_file_path.clone();
                                full_target_file_path.push_str(&target_file_name);

                                if Path::new(&full_target_file_path).exists() {
                                    println!("{} Already Exists, Skipping...", target_file_name);
                                } else {
                                    let data = dataset_rnmd_file.to_str().unwrap();
                                    let mut cmd = Command::new("/bin/bash");
                                    cmd.arg("-c")
                                        .arg(format!(
                                            "python3 {}/bin/gen_targets.py \
                                            {} {} {} {} {} {} {} {} \
                                            {} {} {} {} {} {} {} {} {} \
                                            {} {} {} {} {}",
                                            current_working_dir.clone(),
                                            data.to_string(),
                                            the_ne1_atom.to_string(),
                                            n_atom.to_string(),
                                            internal_atoms.to_string(),
                                            H_BOND_CUT_OFF,
                                            target_file_name,
                                            input_pdb_given.to_string(),
                                            current_fragment,
                                            carved_frag_info[0].to_string(),
                                            carved_frag_info[1].to_string(),
                                            carved_frag_info[2].to_string(),
                                            carved_frag_info[3].to_string(),
                                            carved_frag_info[5].to_string(),
                                            node_id.to_string(),
                                            align_query[0].to_string(),
                                            align_query[6].to_string(),
                                            align_query[4].to_string(),
                                            model_idx,
                                            chain_idx,
                                            residue_idx
                                        ))
                                        .output()
                                        .expect("Failed to execute process");
                                    let output = String::from_utf8_lossy(&cmd.output().stdout);
                                    let error = String::from_utf8_lossy(&cmd.output().stderr);
                                    println!("{}", output);







if carved_frag_info == align_query {
    let old_predicted_receptor_structure = parser_get_structure(&align_query[0], &renamed_pdb[0]);
    let (aa_atom, bb_atom, cc_atom, nx_atom, ee_atom) =
        get_atoms_from_structure(&old_predicted_receptor_structure, &align_query);

    let old_template_peptide_structure = parser_get_structure(&align[8][..4], &dataset_renamed_file[0]);
    let (aa_atom_template, bb_atom_template, cc_atom_template, nx_atom_template, ee_atom_template) =
        get_atoms_from_structure(&old_template_peptide_structure, &align_query);

    let aa_bb_cc_template = vec![aa_atom_template, bb_atom_template, cc_atom_template];
    let nx_ee_template = vec![nx_atom_template, ee_atom_template];

    let aa_bb_cc = vec![aa_atom, bb_atom, cc_atom];
    let nx_ee = vec![nx_atom, ee_atom];

    let rmsd = calculate_rmsd(&aa_bb_cc, &nx_ee, &aa_bb_cc_template, &nx_ee_template);

    let mut new_alignment = align.to_owned();
    new_alignment.push(format!("{:.2}", rmsd));
    new_list_of_unique_alignment.push(new_alignment);
}

if !new_list_of_unique_alignment.is_empty() {
    let best_alignment = new_list_of_unique_alignment.iter()
        .min_by_key(|x| (
            x[x.len() - 2].parse::<f64>().unwrap(),
            -x[x.len() - 1].parse::<f64>().unwrap()
        ))
        .unwrap();
    let best_alignment_query = [        best_alignment[0].to_owned(),
        best_alignment[4][0..1].to_owned(),
        best_alignment[4][1..2].to_owned(),
        best_alignment[4][2..].to_owned(),
        best_alignment[5].to_owned(),
        best_alignment[6].to_owned(),
        best_alignment[8][best_alignment[8].len() - 4..].to_owned()
    ];
    println!("{:?}", new_list_of_unique_alignment);
    if best_alignment[-2].parse::<f64>().unwrap() >= 2.5 {
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
        best_alignment[-2]
    );
    // Reports alignment with least RMSD with "passable" overlap.
    //println!("Best alignment details are - PDB ID = {}, RMSD = {}, SO = {}", best_alignment[8][best_alignment[8].len()-4..], best_alignment[-2], best_alignment[-1]);
    //println!("{:?}", best_alignment);
    //println!("{} {} {} {} {} {} {}", best_alignment[0], best_alignment[4][0..1], best_alignment[4][1..2], best_alignment[4][2..], best_alignment[5], best_alignment[8][best_alignment[8].len()-4..], best_alignment[-2], best_alignment[-1]);



