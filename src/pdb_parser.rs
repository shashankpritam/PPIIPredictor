use std::io::{BufRead, BufReader, Write};
use anyhow::{Result as AnyResult, Error as AnyhowError};
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use std::sync::Arc;
use cgmath::{Point3, EuclideanSpace};



#[derive(Debug, Clone)]
pub struct Structure {
    pub models: Vec<Model>,
    pub chains: Vec<Chain>,
}

impl PartialEq for Structure {
    fn eq(&self, other: &Self) -> bool {
        self.models.len() == other.models.len()
            && self.models.iter().zip(other.models.iter()).all(|(a, b)| a == b)
            && self.chains.len() == other.chains.len()
            && self.chains.iter().zip(other.chains.iter()).all(|(a, b)| a == b)
    }
}

impl Structure {
    pub fn new() -> Self {
        Structure {
            models: Vec::new(),
            chains: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Model {
    pub serial_number: isize,
    pub chains: Vec<Chain>,
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.serial_number == other.serial_number
            && self.chains.iter().zip(other.chains.iter()).all(|(a, b)| a == b)
    }
}

impl Eq for Model {}

#[derive(Debug, Clone)]
pub struct Chain {
    pub id: char,
    pub residues: Vec<Residue>,
}

impl PartialEq for Chain {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
            && self.residues.len() == other.residues.len()
            && self.residues.iter().zip(other.residues.iter()).all(|(a, b)| a == b)
    }
}

#[derive(Debug, Clone)]
pub struct Residue {
    pub name: String,
    pub id: isize,
    pub atoms: Vec<Atom>,
}

impl PartialEq for Residue {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.id == other.id
            && self.atoms.iter().zip(other.atoms.iter()).all(|(a, b)| a == b)
    }
}

#[derive(Debug, Clone, Default)]
pub struct Atom {
    pub serial: isize,
    pub name: String,
    pub alt_loc: char,
    pub res_name: String,
    pub chain_id: char,
    pub res_seq: isize,
    pub icode: char,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub occupancy: f32,
    pub temp_factor: f32,
    pub element: String,
    pub charge: String,
}

impl Atom {
    pub fn get_atom_id(&self) -> isize {
        self.serial
    }
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
            && self.name == other.name
            && self.alt_loc == other.alt_loc
            && self.res_name == other.res_name
            && self.chain_id == other.chain_id
            && self.res_seq == other.res_seq
            && self.icode == other.icode
            && (self.x - other.x).abs() < f32::EPSILON
            && (self.y - other.y).abs() < f32::EPSILON
            && (self.z - other.z).abs() < f32::EPSILON
            && (self.occupancy - other.occupancy).abs() < f32::EPSILON
            && (self.temp_factor - other.temp_factor).abs() < f32::EPSILON
            && self.element == other.element
            && self.charge == other.charge
    }
}

pub fn parse_atom(line: &str) -> Atom {
    let serial: isize = line[6..11].trim().parse().unwrap();
    let name = line[12..16].trim().to_string();
    let alt_loc = line[16..17].chars().next().unwrap_or(' ');
    let res_name = line[17..20].trim().to_string();
    let chain_id = line[21..22].chars().next().unwrap_or(' ');
    let res_seq: isize = line[22..26].trim().parse().unwrap();
    let icode = line[26..27].chars().next().unwrap_or(' ');
    let x: f32 = line[30..38].trim().parse().unwrap();
    let y: f32 = line[38..46].trim().parse().unwrap();
    let z: f32 = line[46..54].trim().parse().unwrap();
    let occupancy: f32 = line[54..60].trim().parse().unwrap_or(0.0);
    let temp_factor: f32 = line[60..66].trim().parse().unwrap_or(0.0);
    let element = line[76..78].trim().to_string();
    let charge = line[78..80].trim().to_string();

    Atom {
        serial,
        name,
        alt_loc,
        res_name,
        chain_id,
        res_seq,
        icode,
        x,
        y,
        z,
        occupancy,
        temp_factor,
        element,
        charge,
    }
}


impl Residue {
    pub fn get_atom_by_name(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|atom| atom.name == name)
    }
}

impl Model {
    pub fn get_atoms(&self) -> Vec<&Atom> {
        self.chains
            .iter()
            .flat_map(|chain| chain.residues.iter())
            .flat_map(|residue| residue.atoms.iter())
            .collect()
    }

    pub fn get_atom_by_name(&self, name: &str) -> Option<&Atom> {
        self.get_atoms().into_iter().find(|&atom| atom.name == name)
    }
}




pub fn parse_pdb_file(file_content: &str) -> AnyResult<Structure> {
    let reader = BufReader::new(file_content.as_bytes());
    let mut structure = Structure::new();
    let mut current_model = Model {
        serial_number: 1,
        chains: Vec::new(),
    };

    let mut current_chain = Chain {
        id: ' ',
        residues: Vec::new(),
    };

    let mut current_residue = Residue {
        name: String::new(),
        id: 0,
        atoms: Vec::new(),
    };

    for result_line in reader.lines() {
        let line = result_line.map_err(AnyhowError::new)?;
        let record_type = line.chars().take(6).collect::<String>();

        match record_type.as_ref() {
            // A new model is encountered
            "MODEL " => {
                current_model = Model {
                    serial_number: line
                        .chars()
                        .skip(10)
                        .take(4)
                        .collect::<String>()
                        .trim()
                        .parse::<isize>()
                        .map_err(AnyhowError::new)?,
                    chains: Vec::new(),
                };
                //println!("New model encountered: {:?}", current_model);
            }
            // An atom is encountered
            "ATOM  " | "HETATM" => {
                let atom = parse_atom(&line);

                // Check if the current residue exists in the current chain
                if let Some(index) = current_chain.residues.iter().position(|r| r.id == current_residue.id) {
                    // If the residue exists, add the atom to it
                    current_chain.residues[index].atoms.push(atom.clone());
                } else {
                    // If the residue does not exist, create a new residue and add the atom to it
                    current_residue.atoms.push(atom.clone());
                    current_chain.residues.push(current_residue.clone());
                }

                // Update the current residue's ID and name
                current_residue.id = atom.res_seq;
                current_residue.name = atom.res_name.clone();
            }
            // A chain terminator is encountered
            "TER   " => {
                current_chain.residues.push(current_residue.clone());
                current_model.chains.push(current_chain.clone());
                //println!("TER encountered: current_chain: {:?}", current_chain);
                current_chain = Chain {
                    id: ' ',
                    residues: Vec::new(),
                };
                //println!("Resetting current_chain: {:?}", current_chain);
            }
            "ENDMDL" => {
                current_chain.residues.push(current_residue.clone());
                current_model.chains.push(current_chain.clone());
                //println!("ENDMDL encountered: current_model: {:?}", current_model);

                // Push the current model to the structure
                structure.models.push(current_model.clone());

                // Reset the current_residue, current_chain, and current_model after pushing them to the structure
                current_residue = Residue {
                    name: String::new(),
                    id: 0,
                    atoms: Vec::new(),
                };
                current_chain = Chain {
                    id: ' ',
                    residues: Vec::new(),
                };
                current_model = Model {
                    serial_number: current_model.serial_number + 1,
                    chains: Vec::new(),
                };
                //println!("Resetting current_residue, current_chain, and current_model");
            }
            _ => (),
        } // <- Add this closing brace
    }

    // Check if the current residue has any atoms and add it to the current chain
    if !current_residue.atoms.is_empty() {
        current_chain.residues.push(current_residue.clone());
    }

    // Check if the current chain has any residues and add it to the current model
    if !current_chain.residues.is_empty() {
        current_model.chains.push(current_chain.clone());
    }

    // Check if the current model has any chains and add it to the structure
    if !current_model.chains.is_empty() {
        structure.models.push(current_model.clone());
    }

    Ok(structure)
}







pub fn write_pdb(structure: &Structure, writer: &mut dyn Write) -> AnyResult<()> {
    for model in &structure.models {
        write!(writer, "MODEL {:>4}\n", model.serial_number)?;
        for chain in &model.chains {
            write!(writer, "CHAIN {:>2} \n", chain.id)?;
            for residue in &chain.residues {
                write!(writer, "RESIDUE {:>3} {:<3} {:>4}\n", residue.id, residue.name, chain.id)?;
                for atom in &residue.atoms {
                    write!(
                        writer,
                        "ATOM  {:>5} {:<4} {:>3} {} {}{:>4}    {:>8.3} {:>8.3} {:>8.3}  {:>6.2} {:>6.2}          {:>2}{:>2}\n",
                        atom.serial,
                        atom.name,
                        residue.name,
                        chain.id,
                        residue.id,
                        atom.icode,
                        atom.x,
                        atom.y,
                        atom.z,
                        atom.occupancy,
                        atom.temp_factor,
                        atom.element,
                        atom.charge
                    )?;
                }
            }
        }
        write!(writer, "ENDMDL\n")?;
    }
    Ok(())
}





pub struct NeighborSearch {
    kdtree: KdTree<f64, Arc<Atom>, [f64; 3]>,
}

impl NeighborSearch {
    pub fn new(atoms: &[Atom]) -> Self {
        let coords: Vec<[f64; 3]> = atoms
            .iter()
            .map(|atom| [atom.x as f64, atom.y as f64, atom.z as f64])
            .collect();
        let atoms: Vec<Arc<Atom>> = atoms.iter().map(|a| Arc::new(a.clone())).collect();
        let mut kdtree = KdTree::new(3);

        for (index, coord) in coords.iter().enumerate() {
            kdtree.add(*coord, atoms[index].clone()).unwrap();
        }

        Self { kdtree }
    }

    pub fn search_neighbors(&self, query: &Atom, radius: f64) -> Vec<Arc<Atom>> {
        let query_coords: [f64; 3] = [query.x as f64, query.y as f64, query.z as f64];
        let neighbor_indices = self.kdtree
            .within(&query_coords, radius.powi(2), &squared_euclidean)
            .unwrap()
            .iter()
            .map(|item| item.1.clone())
            .collect::<Vec<_>>();

        neighbor_indices
    }
}

fn calculate_distance(atom1: &Atom, atom2: &Atom) -> f32 {
    let dx = atom1.x - atom2.x;
    let dy = atom1.y - atom2.y;
    let dz = atom1.z - atom2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn calculate_angle(atom1: &Atom, atom2: &Atom, atom3: &Atom) -> f32 {
    let dx1 = atom1.x - atom2.x;
    let dy1 = atom1.y - atom2.y;
    let dz1 = atom1.z - atom2.z;

    let dx2 = atom3.x - atom2.x;
    let dy2 = atom3.y - atom2.y;
    let dz2 = atom3.z - atom2.z;

    let dot_product = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
    let magnitude1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1).sqrt();
    let magnitude2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2).sqrt();

    let cos_angle = dot_product / (magnitude1 * magnitude2);
    cos_angle.acos().to_degrees()
}
