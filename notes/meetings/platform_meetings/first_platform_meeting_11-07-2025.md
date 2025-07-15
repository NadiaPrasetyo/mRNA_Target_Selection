## Testing Rates and Protocols

- **Testing throughput:**  
    - 1 round per month (or up to 2 per month: 40 candidates tested per month)
    - 10–20 candidates per round can be safely tested using in vivo methods
- **Sequence limits:**  
    - Safe high limit: one sequence (COVID spike, ~4000 kb?)
- **Dosing:**  
    - 2.5 micrograms of one vaccine into mice  
    - Can co-immunize with 2 different constructs at the same time

---

## Team and Affiliations

- **Alex:** Subcontract only until this year for the target selection pillar  
    - AI components seem to be outgoing the pillar—may be of interest for other parts of the platform  
    - Only Holly is working with Alex's team

- **University affiliations:**
    - Wayne Patrick: Victoria University
    - Zohaib, Paul, Nadia: University of Otago
    - Alex Gavryushkin, Holly Pinkney: FoldAI (no university affiliations)
    - Lisa Connor: Malaghan Institute (no university affiliations)
    - Sarah Diermeier: BioTech (no university affiliations)

---

## Projects and Indications

- **Q:** Indications and flagships are still not confirmed for the overall platform?  
    **A:** While there are still decisions to be made across the platform, we can start and learn from small projects to understand the way forward across the platform later on.

---

## Gaps Identified in Platform Research

- Consensus vs. ancestry sequences: flexibility to ask and answer these questions
- Need for consensus sequence per strain vs. building diversity using a mosaic of different strains
- Build a more express high-throughput testing system: serum library, select appropriate sequences, in vitro, mice, many different constructs
- Many vaccines are inactivated viruses (various recurring seasonal strains)—all proteins included
    - Key antibodies: haemagglutinin (for inactivation), but these are not long-lived
    - Goal: engineer sequences (use conserved regions), use different heads for conserved stalks to initiate a longer-term, stronger B cell response against conserved regions while maintaining specificity
- Use AI to predict antigen structures that produce the most stable, surface expression (issue: reliability of these tools)
- Sequence design for efficient protein expression and maximizing product purity for production
- Focus on HA first (influenza)

---

## Goals

1. Engineer B-cell interacting and activating antigen sequences, considering stability, strain, and cross-reactivity against many different strains.
2. Analyze influenza phylogenetic trees to create an ancestral sequence for a vaccine that can protect against various strains, targeting epidemic/pandemic scenarios.

---

## Additional Notes

- There is a group in Auckland active in creating bioinformatics tools to create phylogenetic structures.
- Should we use linear-based and conformational-based predictions?  
    - B-cell conformation is very important; T-cell responses also matter.
    - Conformational epitopes are much more immunogenic.
    - Current focus is on domains and whole antigens, not specific epitopes.
    - Reliability of prediction tools for antigen processing and presentation is a concern.
    - Using only monomers: immunogenicity and flexibility reduce; need to build back T cell epitopes. Narrowing specificity increases, but repertoire is restricted.

---

## Vaccine Design Focus

- Are we building a vaccine against a known or unknown pathogen (influenza specifically)?
    - **A:** The focus is to create a vaccine that can prepare against an unknown strain of a known pathogen, aiming to reduce hospital load and limit spread (full prevention may not be possible).
- Structures must be considered: fix the interfaces but vary the residues (hetero vs. homo trimeric).
- **Current work progress:**  
    - Splitting by domain: head, stalk

