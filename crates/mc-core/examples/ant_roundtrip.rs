//! Writes every SBML model in a folder as Antimony and parses the result
//! back, reporting models whose Antimony has errors.
//!
//! cargo run --release -p mc-core --example ant_roundtrip -- FOLDER

fn main() {
    let dir = std::env::args()
        .nth(1)
        .expect("usage: ant_roundtrip FOLDER");
    let mut files: Vec<_> = std::fs::read_dir(&dir)
        .expect("cannot read folder")
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|e| e == "xml"))
        .collect();
    files.sort();
    let (mut ok, mut bad) = (0, 0);
    for f in &files {
        let Ok(xml) = std::fs::read_to_string(f) else {
            continue;
        };
        let Ok(read) = sbml_io::read(&xml) else {
            continue;
        };
        let Some(model) = &read.document.model else {
            continue;
        };
        let text = mc_core::antimony_out::to_antimony(model, &read.document.namespaces).text;
        let loaded = antimony::parse(&text);
        if loaded.has_errors() {
            bad += 1;
            if bad <= 400 {
                let first = loaded
                    .diagnostics
                    .iter()
                    .find(|d| d.to_string().contains("error"))
                    .map(|d| d.to_string())
                    .unwrap_or_default();
                println!("{}: {first}", f.file_name().unwrap().to_string_lossy());
            }
        } else {
            ok += 1;
        }
    }
    println!("{ok} parse cleanly, {bad} with errors");
}
