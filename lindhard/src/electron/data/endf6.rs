//! A minimal reader for the ENDF-6 records used by File 28.
//!
//! Written from the ENDF-6 Formats Manual (ENDF-102, 2018): an 80-column
//! card has six 11-column data fields (columns 1-66), then MAT (67-70), MF
//! (71-72) and MT (73-75). Floating-point fields may omit the `E` of the
//! exponent (`1.234567+2`); the EPICS files also use Fortran `D` exponents.
//! A blank field is zero.

use super::{AtomBindings, ElectronDataError, Result, ShellBinding, Subshell};

struct Card<'a> {
    line_no: usize,
    fields: &'a str,
    mat: i64,
}

fn parse_err(line_no: usize, reason: impl std::fmt::Display) -> ElectronDataError {
    ElectronDataError::Parse(format!("ENDF line {line_no}: {reason}"))
}

/// An ENDF floating-point field: blank is zero, `D` exponents are
/// accepted, and the `E` of the exponent may be omitted.
pub(super) fn parse_float(raw: &str) -> Option<f64> {
    let raw = raw.trim();
    if raw.is_empty() {
        return Some(0.0);
    }
    let s = raw.replace(['D', 'd'], "E");
    if let Ok(v) = s.parse::<f64>() {
        return Some(v);
    }
    // E-less exponent: the last sign that is not the leading one.
    let p = s[1..].rfind(['+', '-'])? + 1;
    let (m, e) = s.split_at(p);
    format!("{m}E{e}").parse::<f64>().ok()
}

fn float(card: &Card, k: usize) -> Result<f64> {
    let raw = &card.fields[11 * k..11 * (k + 1)];
    parse_float(raw).ok_or_else(|| {
        parse_err(
            card.line_no,
            format!("field {} is not a number: {:?}", k + 1, raw.trim()),
        )
    })
}

fn int(card: &Card, k: usize) -> Result<i64> {
    let raw = card.fields[11 * k..11 * (k + 1)].trim();
    if raw.is_empty() {
        return Ok(0);
    }
    raw.parse::<i64>().map_err(|_| {
        parse_err(
            card.line_no,
            format!("field {} is not an integer: {raw:?}", k + 1),
        )
    })
}

/// The MF = 28, MT = 533 cards of a file, in order.
fn mf28_cards(text: &str) -> Result<Vec<Card<'_>>> {
    let mut cards = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = line.trim_end_matches('\r');
        if line.len() < 75 {
            continue; // too short to carry MF/MT: not a File 28 card
        }
        if !line.is_ascii() {
            return Err(parse_err(line_no, "non-ASCII characters"));
        }
        let (mf, mt) = (line[70..72].trim(), line[72..75].trim());
        if mf != "28" || mt != "533" {
            continue;
        }
        let mat = line[66..70]
            .trim()
            .parse::<i64>()
            .map_err(|_| parse_err(line_no, "MAT field is not an integer"))?;
        cards.push(Card {
            line_no,
            fields: &line[..66],
            mat,
        });
    }
    Ok(cards)
}

pub(super) fn read_mf28(text: &str) -> Result<Vec<AtomBindings>> {
    let cards = mf28_cards(text)?;
    if cards.is_empty() {
        return Err(ElectronDataError::Parse(
            "no ENDF-6 MF=28 MT=533 records found".into(),
        ));
    }
    let mut atoms = Vec::new();
    let mut start = 0;
    while start < cards.len() {
        let mat = cards[start].mat;
        let end = start
            + cards[start..]
                .iter()
                .position(|c| c.mat != mat)
                .unwrap_or(cards.len() - start);
        if let Some(atom) = read_section(&cards[start..end])? {
            atoms.push(atom);
        }
        start = end;
    }
    Ok(atoms)
}

/// One material's section: HEAD then NSS subshell LIST records.
fn read_section(cards: &[Card]) -> Result<Option<AtomBindings>> {
    let head = &cards[0];
    let za = float(head, 0)?;
    if !(za.is_finite() && za >= 1000.0 && za.fract() == 0.0 && za % 1000.0 == 0.0) {
        return Err(parse_err(
            head.line_no,
            format!("ZA = {za} is not an elemental ZA (1000 Z)"),
        ));
    }
    let z = za / 1000.0;
    let nss = int(head, 4)?;
    if nss < 1 {
        return Err(parse_err(head.line_no, format!("NSS = {nss}, need >= 1")));
    }
    let mut next = 1;
    let mut shells = Vec::new();
    for _ in 0..nss {
        let Some(h) = cards.get(next) else {
            return Err(parse_err(
                cards[cards.len() - 1].line_no,
                "section ends before all NSS subshells were read",
            ));
        };
        let subi = float(h, 0)?;
        let nw = int(h, 4)?;
        let ntr = int(h, 5)?;
        if ntr < 0 || nw != 6 * (1 + ntr) {
            return Err(parse_err(
                h.line_no,
                format!("NW = {nw} with NTR = {ntr}; expected NW = 6 (1 + NTR)"),
            ));
        }
        if !(subi.fract() == 0.0 && (1.0..=f64::from(Subshell::MAX_DESIGNATOR)).contains(&subi)) {
            return Err(parse_err(
                h.line_no,
                format!("SUBI = {subi} is not a subshell designator"),
            ));
        }
        let data_lines = usize::try_from(nw / 6).unwrap_or(0);
        let Some(first) = cards.get(next + 1) else {
            return Err(parse_err(h.line_no, "LIST record has no data line"));
        };
        if cards.len() < next + 1 + data_lines {
            return Err(parse_err(h.line_no, "LIST record is truncated"));
        }
        let (ebi, eln) = (float(first, 0)?, float(first, 1)?);
        // Validate the skipped transition lines as numbers too, so a
        // corrupt record is not silently accepted.
        for c in &cards[next + 1..next + 1 + data_lines] {
            for k in 0..6 {
                float(c, k)?;
            }
        }
        next += 1 + data_lines;
        shells.push((subi as u8, ebi, eln, h.line_no));
    }
    if next != cards.len() {
        return Err(parse_err(
            cards[next].line_no,
            "unexpected records after the last subshell",
        ));
    }
    if z > crate::elements::NUM_ELEMENTS as f64 {
        return Ok(None);
    }
    let shells = shells
        .into_iter()
        .map(|(d, ebi, eln, line_no)| {
            ShellBinding::new(Subshell::new(d)?, ebi, eln).map_err(|e| parse_err(line_no, e))
        })
        .collect::<Result<Vec<_>>>()?;
    AtomBindings::new(z as u8, shells)
        .map(Some)
        .map_err(|e| parse_err(head.line_no, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn endf_floats_with_and_without_exponent_letter() {
        let cases = [
            (" 1.500000+3", Some(1500.0)),
            ("1.0000D+11", Some(1.0e11)),
            ("13.6000000", Some(13.6)),
            (" .999241400", Some(0.9992414)),
            ("-2.5-1", Some(-0.25)),
            ("", Some(0.0)),
            ("7", Some(7.0)),
            ("1.2.3", None),
        ];
        for (raw, want) in cases {
            assert_eq!(parse_float(raw), want, "{raw:?}");
        }
    }
}
