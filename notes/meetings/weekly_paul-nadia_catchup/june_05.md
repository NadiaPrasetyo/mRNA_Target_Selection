To ask:
- Are we using this academically only? some of the tools allow downloading only for academic purposes not comercial:
  - For now (non commercial)
- If the goal is to create mRNA vaccines, do the antigens and epitopes we find have to be novel and unpattented?
  - No for now, if you find something we then check


To report:
- Zohaib is working on an experimental pipeline for the universal flu and extending to s. aureus atm
- I am working on using epitope predictive tools to figure out good epitopes from the antigens - which can be used for training data later
- I will continue with getting the metrics of a good target for vaccine - antigenicity, allergenicity, immunogenicity, population coverage, and conservation
- Looking throught the antigens we fetched, some of them also had domains/features/function annotations that could be looked at

Allen G. Rodrigo - great paper on why a consensus sequence is a terrible idea for influenza vaccine
    -> consensus sequence = non biological sequence that maybe doesn't even fold - unlikely to be a very effective vaccine
    -> viral evolutionary trees 

Do instead ancestral sequence reconstruction - run evoution back in time: real sequences (or closer to reality)

quality metrics for mapping: pident, bit score, 

alignment: calculate conservation - tools available

Other opportunities:
**Biological White-Hat Project**
The Association of Biosafety in Australia and New Zealand: can we design a better pathogen?
-> can we detect if anyone using AI for malicious use to develop sequences e.g for pathogens
old school methods: homology mapping
using hidden markov models (mmseq2, hammer) to generate sequences
biocontrol: NZ has a long history with using animals, pathogens to control 
maximatosis, kalesii virus -> control rabbit population
-> designing pathogens

rather than limiting- > better detection: building database of harmful biological sequences - considering how quickly things evolve, how would you keep it updated:
    - traditionally: hire curaters that check the literature adn data surrounding the sequences
    - automated using AI: transformer models - check the relation between the an unknown sequence and a database
  
Essential penetration testing
Help paul be a white hat hacker

PRO: public research organization