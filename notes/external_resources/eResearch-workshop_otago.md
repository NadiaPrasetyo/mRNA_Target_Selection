# eResearch @ Otago

## Updates on NZ

- NeSI has been absorbed into REANNZ (merged national scale facilities)
- AI intelligence has a seen effect on HPC infrastructure, open source capabilities, GPU demand
- Otago compute centralization has boosted mid-tier compute
- RTIS facilities and support has grown

**Website:** [https://www.otago.ac.nz/eresearch/](https://www.otago.ac.nz/eresearch/)

---

## Speaks

### eResearch facilities and services @ Otago - RTIS support

Research focused support but has + running old systems that are slowly migrated to the GPU and centralized systems; around 20 people team that takes care of equipment, solutions, and programming

#### Major systems in Otago

- **High Capacity Storage (HCS):** primary location to store data on campus
- **HCS - S3 capabilities:** cloud work, local based, network connectivity (about to be launched)
- **Research Powerscale:** (higher connection and speed) 3 petabyte
- **Weka:** small chunk HPC scale store systems - 4 Tb

#### Data transfer mechanism

- **Globus:** (encrypted, high speed, higher technical requirement): best for big data and out of campus
- **Otago File-Transfer:** smaller external stuff (no need permissions, just login)
- **REANNZ FileSender service:** (allow sending 5TB of data)

#### Datacentre Infrastructure at Otago

- 43 Node clusters: large CPU, large memory with around 50 GPU, 5000+ cores
- lots of apps associated - on the web, not in desktop (can request specific apps)

#### Support & Training

- collaborate with ResBaz, the national carpentries
- can create/tailor training as required
- carpentries course
- application installation, etc

#### Consultancy

- get stuff and systems set up
- data management,
- big data and big data consideration

---

### Outside of Otago scope of services: REANNZ

NeSI is a collaboration around high performance computing - has now merged into REANNZ (Crown owned entity - provide research and network capabilities)

> we all use REANNZ via the internet that we use in the university for research and develop services

**NeSI services within REANNZ:**

- access to HPC supercomputers, large clusters
- training
- genoa nodes, milan nodes, kit running
- GPUs (A100, H100s, scale up)
- support team: based on campus and around the country

**Involved in community building aspects:**

- national training
- research software engineering community
- career progression advise

---

## Surfing the Wave: What is AI’s Role in Qualitative Health Research?

- Spirituality and AI: spiritual care in NZ in hospitals and clinical care (religious, secular, existential)
- AI has come along as a useful tool for qualitative studies - reduce labout surrounding transcription, summaries and coding, literature reviews
- includes: Whisper X, MAXQDA, Copilot, Perplexity, NotebookLM, Gemini, etc
- AI used for faster transcription and analysis, support for writing and editing, enhanced thematic exploration

**Ethical challenges:** data privacy, confidentiality, data sovereignty (especially for Maori data), bias and inequity in AI models, transparency issues

**Environmental and cultural considerations:** many questions and discussions surrounding energy consumption, te tiriti

**Q:** how to use copilot, gemini for literature review?  
**A:** I am learning how to use the prompts and questions to get the results that is good: prompt engineering - can also use AI tools to make good AI prompts. AI can help refine ideas and writings that we as researchers do a lot of.

**Q:** For thematic coding, does the researcher read every single results (line by line) or just use with AI.  
**A:** It's interesting how it uses language but so far I have only done coding by hand as a human - and refine or tesitng using AI.

---

## Leveraging Gen-AI for teaching and learning

**NEST: AI agents platform Management: 4 parts**

1. NEST-DB
2. NEST-Gateway - acts as an API (open AI compatible endpoint)
3. NEST-Auth
4. NEST-Admin

**Future:**

5. NEST-Chat
6. NEST-Embed
7. NEST-MCP
8. NEST-RAG (mix of vector and graph databases)

- divide users into roles permissions and groups
- models have different providers (including locally hosted), different access permissions (data confidentiality), creating agents specific to its use (easily accessible, used)
- used across the university for teaching, learning, supports

**Q:** is NEST already available for staff?  
**A:** it's in limited trial at a time - helping users with agent building, prompt engineering, for educational outcomes

**Q:** What are your use cases at the moment?  
**A:** teaching and education uses: constructive feedbacks, learning skills feedback, planned to deploy for second semester with some courses. e.g providing feedback on writing levels.

**Q:** does it allow audio input?  
**A:** no, we are optimizing user experience; we are waiting for better audio input systems to be developed.

**Q:** Does NEST keep data in Otago systems/ local by choosing the University of Otago provider?  
**A:** Yes, it'll maintain the data sovereignty.

---

## Complex history of invasive Rattus in New Zealand (Aotearoa) revealed by whole genomic resequencing

Mainly use NeSI and Aoraki clusters to analyze the population data of invasive black rats and european brown rats in New Zealand.

**Management methods in NZ:**

- Physical snap traps - difficult to scale up, relatively humane
- Poisoning - off target impacts
- Genetic (Gene drive) - species specific, but potential for resistance and spread introduction

**Gene drive technology:** use T haplotype using super mendelian inheritance - spread to 90% offsprings (preferential inheritance), linked to ingertility genes (Tiam2 & SMOK)

20 year time span for trapping and collection of samples. Whole genome sequences: procide more insight with 7 singletons observed across NZ

> use modeling to see and predict how the gene drive is going to affect the different population clusters of rats across nz.

---

## Modelling ozone concentration with Google’s GraphCast deep learning model

- developing deep learning models for various things - e.g. modeling ozone levels and how well the ozone hole is recovering. Findings: it's not recovering as well as we thought and the hole is bigger than it should be.
- next steps: why? - solution: build a model that can analyze and model ozone levels across the world to understand what causes the holes, predict what's happening based on daily measurements via satellites across the globe. Include data with water vapor, atmosphere composition, and water temperature.

**Base model:** GraphCast - model that successfully model weather with faster and more accurate prediction. Model is based on an icosahedral mesh around the globe with edges and nodes with each has a tiny neural network that passes the messages (next level convolutional network): number of times, mimic the swirling of the gasses and maps back to the mesh; take into account the column of different gases to run prediction a few days forward to see what is happening.

**Compute used:** Aoraki clusters:  
1 GPU (H100) - 72.1 GPU RAM  
trained on 2004-2017 data, valiidated on 2018-2020 data  
~5hr per 1 training epoch (93epoch)

**results:** ozone predicted on different longiture and latitudes

**Q:** Is the chemicals a discrete list that you are looking for to see their impacts on the ozone? will there be a potential missed confounding chemical that is not detected?
**A:** Yes it is a discrete list at this point and it is a limiation, especially with missing factors that may be critical.

---

## Unlocking Research Potential: Leveraging Dimensions Analytics to build new collaborations

**Dimensions:** access via the library - see what is happening in a particular research discipline  
Platform that connects all information regarding a particular researcher, research field, datasets, publications, patents, grants, collaborative work, research topics, etc

- can see what opportunities are available in sustainability goals, field of research, keywords/concepts, where research is happening
- commercializing research: see what organizations are working with research organizations, what patents are being made, what organizations are assigning research.
- can see analytical views on publications: citations over time, publications over time, FCR (field normalized citation metric, etc)

**Q:** How to access it from otago website?  
**A:** Under library databases: Dimensions (digital science) - log in using otago credentials

**Q:** How often do you update information and where from?  
**A:** Crossref, pubmed, 13 direct connections with organizations - depends how quick it is updated to dimensions but 2-3 days for crossref. Funding information directly from the funding provider, so no latency from 3-5 years (no need to be confirmed via publications)

**Q:** Can it be used to find a potential supervisor for projects/postgrad research?  
**A:** What I would recomment: use granular search, FOR, networks, and see how established they are. From there what you can actually get emails/google search, see if they are taking opportunities for projects/research

---
## Managing Research Data with the New Library RDM Module and Updates on Otago Research Information Systems

**New Research Data Management Module (RDM):**  
Practical guideline - mainly designed for post-graduate students and researchers  
→ Created to improve awareness of good RDM practices among project students and new researchers  
→ 4 parts, 12 lessons and include key questions, checklists, tips, visualisations, real examples and explanation, and external resources

[Module link](https://otago.libguides.com/data_management)

**Otago Research Information Systems: Conduit**  
Prepare, submit human ethics submissions platform. It also does all the assessments and runs meetings in the administration of human ethics application.

**Otago Research Information Systems: Research Administration System (RAS):**  
Administration platform accessed by all academics, based on their roles have different levels of access. Includes budget builder capabilities.

**In development: Ngai Tahu Research Consultation Committee**  
Process for consultation with Māori, replace the NTRCC process via the ORIS platform. Includes the application, submission, administration and assessments, and communication.

**In development: Research Committee Funded Projects**  
Process for applying for Research Committee funding, includes application, assessment, notification processes via RAS. Will be used for various grants.

---

## Aotearoa Genomic Data Repository (AGDR): A Resource for Researchers and Kaitiaki

**Genomics Aotearoa:**  
MBIE funded platform to support developments in genomics and bioinformatics - received funding to continue 2025-2030 (in transition stage).

**AGDR:**  
Has been co-developed by Genomics Aotearoa and NeSI - has been running for 5 years.  
Main goal: a national genomics data repository including bespoke resources.....(FINISH LATER). Focus on Taonga species (indigenous and endemic flora and fauna).

Data management requirement of research data: **FAIR**
- **F**indable
- **A**ccessible
- **I**nteroperable
- **R**eusable

For indigenous data: **CARE**
- **C**ollective benefit
- **A**uthority to control
- **R**esponsibility
- **E**thics

Contains an AGDR Māori consultation board (Kaitiaki consultation for approval for all genomics projects) to help with research involving indigenous data and complying with the Te Tiriti policies.

---

## Adolescence After the Algorithmic Turn: The Autonomy-Enhancing and Autonomy-Inhibiting Impacts of Social Media

Impacts of technology on adolescent development.

**Glossary:**
- **Affective experiences:** recognition, misrecognition, and non-recognition (from others)
- **Hybrid social worlds:** both online and offline identity development; online are not less real than offline world
- **Political identities:** based on cancel culture, etc

**Factors findings:**
- **Internal factors:** controlled interactions between self and social media, social network, expectations, and communities
- **External factors:** out of control (platform affordances, sorting logics, filtration based on interest) - limits the potential identities available to certain people

**Q:** What was the response towards the governmental social media bans?  
**A:** There is interesting point that was not discussed in the study but an observation in a following studies: Young people are shutting down and deleting social media themselves.

**Q:** Do the young people realize the paths of social media engagement?  
**A:** Yes, they realize the algorithms and bots in social media is pushing certain content towards specific users; they are attuned to understanding the social media better than some of the naive adults.

---

## Cybersecurity for Research and Researchers at Otago

**Threat Landscape for Cybersecurity in NZ and Globally:**
- NCSC annual report: increasing number of cyber threat sources
- CyberCX annual report: deteriorating cyber threat landscape, increased cyber attacks
- Otago is NOT immune - we had more than a dozen cybersecurity incidents and thousands of events, breaches that reduce the services or data lost

**Phishing M365 Attacks:**  
Target the Microsoft 365 accounts, which is the basis of Otago resources

**Threat landscape for research:**
- Hidden foreign military purpose
- Military & Intel risk
- Sovereign data risk
- Reputational risk
- Non-civilian research landscape/proposal

**Access permissions classifications:**
1. **Public** - public release and domain  
    → e.g. Public published research
2. **Internal Use** - university operations but no personal or sensitive business information
3. **Business Confidential** - non-public strategic/commercial business information  
    → e.g. Microsoft Teams (Otago Uni's secured instance of Teams)
4. **Private Confidential** - non-public personal information about any identifiable individuals  
    → HCS, Virtual machines, aggregated and all PI removed before collaboration phase, hospital data
5. **Restricted** - very sensitive business or personal information

**Top tips:**
- Noone's email
- Training in Blackboard: content/staff cyber security training
- Classify your data and follow guidance
- Report incidents or issues to cybersecurity@otago.ac.nz

---

## The Impact of New Zealand's Changing Climate on Agriculture and Renewable Energy Sectors

**Active projects:**  
To model climate projections - link climate variables to environmental data to understand the climate impacts including drought on agriculture and renewable energy.

**Modeling of climate change:**  
Low resolution columns on the atmosphere and in the ocean - get a low resolution global climate landscape + regional climate model: higher resolution data on specific region.

**NEW MODEL:**  
Conformal cubic atmospheric model developed by CSIRO. High resolution but also high compute.

Downscaling GCMs for climate impact studies to reduce the time intensity and resource intensive: approaches include using regional climate models, and machine learning regressions - generate large ensemble that is a better representation of the data.

---

## Assembling Reference Genomes to Uncover Population Histories of Endangered Dolphins

Genome repository for endangered dolphins: Hector's and Maui dolphins - endemic to New Zealand, amongst the smallest of the dolphins. Affected by fisheries impacts.

**Issues:**  
Non-model species, endangered, marine organisms, non-optimal input data → difficult to establish a reference genome. Now published for conservation efforts.

Look into the demographic history, heterozygosity → changes in the population size over time, divergence of Maui and Hector's dolphins happened due to land breach stopping gene flow. Using the populational data and genome wide heterozygosity provides insights to risks to inbreeding and depression in gene diversity.

**The value of these genomic resources for the dolphins:**  
Reference chromosome level genomes with pipeline applicable for other groups. Determined current status as subspecies and lack of genetic flow between Hector's and Maui dolphins + identify population of the Otago dolphins as part of the south population and not the east coast population.

**Q:** Are there any attempt in breeding the south population with other populations to prevent inbreeding and improve the population in such small population?  
**A:** Manipulating dolphins are very hard, and they live in very small clusters. Plus we would need ethical considerations to do the genetic rescue.

---
