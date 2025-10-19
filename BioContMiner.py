import os
from Bio import Entrez, Medline
import pandas as pd
from urllib.error import HTTPError

Entrez.email = "sjdhhe@gmail.com"

def search_pubmed(query, max_results=1000):
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
            print(f"HTTPError: {e} - Skipping  IDs: {chunk_ids}")
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
        country = record.get("PL", "")
        authors = record.get("FAU", [])
        year = record.get("Year", "")
        affiliations = record.get("AD", [])
       
        email_author_pairs = []
        if isinstance(affiliations, list) and isinstance(authors, list):
            for i, author in enumerate(authors):
                author_emails = []
                if i < len(affiliations): 
                    for line in affiliations[i].split(";"):
                        if "@" in line:
                            author_emails.append(line.strip())
                
                if not author_emails:
                    author_emails.append("No emails found")
                
                email_author_pairs.append({"Author": author, "Emails": author_emails})
        
        data.append({
            "PubMed ID": pubmed_id,
            "PubMed Link": pubmed_link,
            "Year": year,
            "Country": country,
            "Author-Email Pairs": email_author_pairs
        })
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
        country = entry["Country"]
        for pair in entry["Author-Email Pairs"]:
            flat_data.append({
                "PubMed ID": pubmed_id,
                "PubMed Link": pubmed_link,
                "Year": year,
                "Country": country,
                "Author": pair["Author"],
                "Emails": ", ".join(pair["Emails"])
            })
  
    df = pd.DataFrame(flat_data)
    save_name = input("Enter the name to save the file")
    df.to_csv(f"{save_name}.csv", index=False)
    print(f"Data saved to {save_name}.csv")

if __name__ == "__main__":
    main()