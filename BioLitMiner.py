import os
from Bio import Entrez, Medline
import pandas as pd
from urllib.error import HTTPError
import re 
from difflib import get_close_matches

Entrez.email = "tutue@gmail.com"

disease_names = [
    'cancer', 'melanoma', 'leukemia', 'lymphoma', 'carcinoma', 'sarcoma', 'brain tumor', 'breast cancer', 'prostate cancer', 'lung cancer', 'colorectal cancer', 'pancreatic cancer',
    'liver cancer', 'ovarian cancer', 'cervical cancer', 'bladder cancer', 'kidney cancer', 'hodgkin lymphoma', 'non-hodgkin lymphoma', 'adenocarcinoma', 'squamous cell carcinoma',
    'basal cell carcinoma', 'osteosarcoma', 'liposarcoma', 'glioblastoma', 'meningioma', 'small cell lung cancer', 'non-small cell lung cancer', 'hepatocellular carcinoma',
    'renal cell carcinoma', 'testicular cancer', 'thyroid cancer', 'esophageal cancer', 'gastric cancer', 'head and neck cancer', 'pancreatic adenocarcinoma', 'multiple myeloma', 'chronic lymphocytic leukemia',
    'diabetes', 'type 1 diabetes', 'type 2 diabetes', 'gestational diabetes', 'prediabetes', 'latent autoimmune diabetes in adults', 'maturity onset diabetes of the young',
    'hypertension', 'primary hypertension', 'secondary hypertension', 'prehypertension', 'white-coat hypertension', 'malignant hypertension','asthma', 'allergic asthma', 'exercise-induced asthma', 'occupational asthma',
    'nocturnal asthma', 'cough-variant asthma', 'brittle asthma', 'arthritis', 'osteoarthritis', 'rheumatoid arthritis', 'psoriatic arthritis', 'gout', 'ankylosing spondylitis', 'juvenile idiopathic arthritis', 'septic arthritis',
    'cardiovascular disease', 'coronary artery disease', 'heart attack', 'angina', 'stroke', 'ischemic stroke', 'hemorrhagic stroke', 'heart failure', 'arrhythmia', 'atrial fibrillation',
    'ventricular tachycardia', 'peripheral artery disease', 'aortic aneurysm', 'cardiomyopathy', 'congenital heart disease', "alzheimer's disease", "early-onset alzheimer's", "late-onset alzheimer's",
    'vascular dementia', 'lewy body dementia', 'frontotemporal dementia', "parkinson's disease", "idiopathic parkinson's disease", 'atypical parkinsonism',
    'multiple system atrophy', 'progressive supranuclear palsy', "early-onset parkinson's", "parkinson's disease dementia",
    'chronic obstructive pulmonary disease', 'chronic bronchitis', 'emphysema', 'refractory asthma', 'kidney disease', 'chronic kidney disease', 'acute kidney injury', 'nephrotic syndrome',
    'polycystic kidney disease', 'glomerulonephritis', 'liver disease', 'cirrhosis', 'hepatitis', 'hepatitis a', 'hepatitis b', 'hepatitis c', 'hepatitis d', 'hepatitis e', 'non-alcoholic fatty liver disease', 'alcoholic liver disease',
    'liver cancer', 'depression', 'major depressive disorder', 'persistent depressive disorder', 'anxiety disorder', 'generalized anxiety disorder', 'panic disorder', 'social anxiety disorder', 'bipolar disorder',
    'schizophrenia', 'obsessive-compulsive disorder', 'post-traumatic stress disorder', 'autoimmune disease', 'lupus', 'multiple sclerosis', "crohn's disease", 'ulcerative colitis',
    'celiac disease', 'infectious disease', 'tuberculosis', 'hiv/aids', 'malaria', 'influenza', 'covid-19', 'pneumonia', 'neurological disorder', 'epilepsy', 'migraine', 'amyotrophic lateral sclerosis', "huntington's disease"
]

organism_terms = [
    'homo sapiens', 'human', 'plant', 'mus musculus', 'mouse', 'rattus norvegicus', 'rat', 'escherichia coli', 'e. coli', 'saccharomyces cerevisiae', 'yeast', 'drosophila melanogaster', 'fruit fly', 'caenorhabditis elegans', 'worm',
    'mycobacterium tuberculosis', 'staphylococcus aureus', 'influenza a virus', 'covid', 'sars-cov-2', 'plasmodium falciparum', 'hepatitis c virus', 'human immunodeficiency virus', 'hiv', 'zika virus', 'ebola virus', 'arabidopsis thaliana', 
    'thale cress','danio rerio', 'zebrafish', 'candida albicans', 'schizosaccharomyces pombe', 'fission yeast', 'neurospora crassa', 'streptococcus pneumoniae', 'pneumococcus', 'salmonella enterica', 'salmonella', 'clostridium difficile', 
    'chlamydia trachomatis', 'aspergillus fumigatus', 'oryza sativa', 'rice', 'zea mays', 'maize', 'corn', 'gallus gallus', 'chicken', 'canis lupus familiaris', 'dog', 'felis catus', 'cat', 'macaca mulatta', 'rhesus macaque', 'pan troglodytes', 'chimpanzee', "Dengue virus", 
    'Zika virus', 'Chikungunya', 'Yellow fever virus', 'Japanese encephalitis virus', 'West Nile virus'
]

sample_terms = [
    'blood', 'urine', 'plasma', 'serum', 'swab', 'tissue', 'biopsy', 'ffpe', 'formalin-fixed paraffin-embedded', 'saliva', 'sputum', 'cerebrospinal fluid', 'csf', 'stool', 'fecal sample', 'synovial fluid', 'bone marrow', 'amniotic fluid', 
    'bronchoalveolar lavage', 'bal', 'cell culture', 'cell lysate', 'dna', 'rna', 'protein extract', 'peripheral blood mononuclear cells', 'pbmc', 'hair', 'nail clippings', 'fresh frozen','skin scrape', 'buccal swab', 'fixed tissue',
    'ascites fluid', 'pleural fluid', 'vitreous humor', 'aqueous humor', 'frozen tissue', 'cryopreserved sample', 'microbiome sample', 'metabolite extract', 'breast milk', 'colostrum', 'fast frozen', 'semen', 'seminal fluid', 'sweat',
    'gastric fluid', 'bile', 'pancreatic fluid', 'lymph node aspirate', 'fine needle aspirate', 'fna', 'organoid', 'spheroid', 'exosome', 'extracellular vesicles', 'tumor microenvironment sample', 'single-cell suspension',
    'paraffin section', 'histological section', 'environmental swab', 'surface swab', 'soil sample', 'water sample', 'plant tissue', 'leaf extract', 'biofilm', 'microbial culture'
]

def search_pubmed(query, max_results=10):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_details(id_list, chunk_size=10):
    records = ""
    for i in range(0, len(id_list), chunk_size):
        chunk_ids = id_list[i:i + chunk_size]
        ids = ",".join(chunk_ids)
        try:
            handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
            chunk_records = handle.read()
            handle.close()
            records += chunk_records
        except HTTPError as e:
            print(f"HTTPError: {e} - Skipping IDs: {chunk_ids}")
            continue
    return records

def parse_records(records):
    handle = records.splitlines(True)
    parsed_records = list(Medline.parse(handle))
    return parsed_records

def extract_info(parsed_records):
    data = []
    for record in parsed_records:
        pubmed_id = record.get("PMID", "")
        pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
        dp_field = record.get("DP", "")
        year_match = re.search(r"\b(19|20)\d{2}\b", dp_field)
        year = year_match.group(0) if year_match else ""
        funding = record.get("GR", [])
        funding_str = ", ".join(funding) if funding else "No funding information"
        
        authors = record.get("FAU", [])
        affiliations = record.get("AD", [])
        
        author_aff_pairs = []
        for i, author in enumerate(authors):
            if i < len(affiliations):
                affiliation = affiliations[i]
            else:
                affiliation = affiliations[-1] if affiliations else ""
            author_aff_pairs.append({"Author": author, "Affiliation": affiliation})
        
        data.append({
            "PubMed ID": pubmed_id,
            "PubMed Link": pubmed_link,
            "Year": year,
            "Author-Aff Pairs": author_aff_pairs,
            "Funding": funding_str
        })
    return data

def extract_aff(text):
    email_pattern = r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}'
    text = re.sub(email_pattern, '', text)
    text = re.sub(r'Electronic address:', '', text)
    parts = [part.strip() for part in text.split(",") if part.strip()]
    
    if len(parts) >= 3:
        institute = parts[0]
        state = parts[-2]
        country = parts[-1]
    elif len(parts) == 2:
        institute = parts[0]
        state = ""
        country = parts[1]
    elif len(parts) == 1:
        institute = parts[0]
        state = ""
        country = ""
    else:
        institute = ""
        state = ""
        country = ""
    return pd.Series([institute, state, country])

def clean(data):
    data["Email Only"] = data["Affiliation"].str.extract(r'([A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,})')
    data[["Institute", "State", "Country"]] = data["Affiliation"].apply(extract_aff)
    data["Email Only"] = data["Email Only"].fillna("No email found")
    return data

def match_terms(text, terms):
    found = set()
    for word in terms:
        if word.lower() in text:
            found.add(word)
    return list(found)

def add_ass_data(data):
    for idx, row in data.iterrows():
        try:
            fetch = Entrez.efetch(db="pubmed", id=row["PubMed ID"], rettype="medline", retmode="text")
            record = list(Medline.parse(fetch))[0]
            abstract = record.get("AB", "").lower()
            data.at[idx, "Diseases"] = ", ".join(match_terms(abstract, disease_names))
            data.at[idx, "Organisms"] = ", ".join(match_terms(abstract, organism_terms))
            data.at[idx, "Samples"] = ", ".join(match_terms(abstract, sample_terms))
        except Exception as e:
            data.at[idx, "Diseases"] = ""
            data.at[idx, "Organisms"] = ""
            data.at[idx, "Samples"] = ""
    return data

def main():
    query = input("Input the query (Use only keywords): ")
    id_list = search_pubmed(query)
    records = fetch_details(id_list)
    parsed_records = parse_records(records)
    data = extract_info(parsed_records)

    flat_data = []
    for entry in data:
        pubmed_id = entry["PubMed ID"]
        pubmed_link = entry["PubMed Link"]
        year = entry["Year"]
        funding = entry["Funding"]
        for pair in entry["Author-Aff Pairs"]:
            flat_data.append({
                "PubMed ID": pubmed_id,
                "PubMed Link": pubmed_link,
                "Year": year,
                "Author": pair["Author"],
                "Affiliation": pair["Affiliation"],
                "Funding": funding
            })

    df = pd.DataFrame(flat_data)
    cl_data = clean(df)
    fin_data = add_ass_data(cl_data)
    fin_data.to_csv("C:/Users/aksha/Desktop/Ext_Data.csv", index=False)
    return fin_data

if __name__ == "__main__":
    main()