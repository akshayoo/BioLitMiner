Biomedical Literature Mining Tools

This repository contains two Python scripts for automated retrieval, processing, and analysis of biomedical literature from PubMed. These tools leverage Biopython and pandas to extract article metadata, authorship information, affiliations, and biomedical annotations.

---------------------------------------------------------

1. BioLitMiner.py

BioLitMiner.py is designed to retrieve and annotate biomedical literature from PubMed based on user-provided keywords. It automates the extraction of detailed article metadata and provides a structured dataset for downstream analysis.

Features

Searches PubMed for articles matching user queries.
Retrieves records in batches to handle large datasets.
Extracts metadata including:
PubMed ID
Publication year
Authors
Affiliations (split into institute, state, country)
Funding information
Annotates abstracts for mentions of:
Diseases
Organisms
Sample types
Flattens author-affiliation pairs into a structured pandas DataFrame.
Handles errors during fetching to ensure continuity.
Exports cleaned and annotated data to a CSV file on the desktop.

Use Cases

Literature surveys
Research tracking
Biomedical data analysis

---------------------------------------------------------

2. BioContMiner.py

BioContMiner.py is a PubMed email extractor that automates the collection of author contact information from scientific publications.

Features

Searches PubMed for articles matching user-defined keywords (up to 1000 results by default).
Fetches articles in chunks to handle large datasets efficiently.
Extracts metadata including:
PubMed ID
Publication year
Country of publication
Authors and their email addresses
Matches authors to their affiliation lines to find emails, defaults to "No emails found" if unavailable.
Flattens the data into a pandas DataFrame, with one row per author.
Exports the dataset to a user-specified CSV file.

Use Cases

Academic networking
Outreach to authors
Building contact lists from scientific publications

---------------------------------------------------------

Installation

Clone the repository:
git clone https://github.com/yourusername/biomedical-literature-miners.git
cd biomedical-literature-miners


Install dependencies:
pip install -r requirements.txt

Dependencies include: Biopython, pandas.

---------------------------------------------------------

Usage

BioLitMiner.py
python BioLitMiner.py

Enter your PubMed search query when prompted.
The script will fetch, annotate, and save a CSV file with article metadata and biomedical annotations.

BioContMiner.py
python BioContMiner.py

Enter your PubMed search query and desired CSV filename.
The script will fetch author metadata and emails, and save it as a CSV.

---------------------------------------------------------

License
This repository is licensed under the MIT License.
