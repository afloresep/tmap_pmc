import os
import pathlib
import re
import xml.etree.ElementTree as ET
from typing import Tuple


def process_article(path: str) -> Tuple:
    with open(path, "r") as f:
        root = ET.parse(f)
        body = root.find("body")
        text = ET.tostring(body, encoding="utf-8", method="text").decode("utf-8")
        text = re.sub("[^a-zA-Z - — —]", " ", text)
        text = re.sub("[— — (\s\s+)]", " ", text)
        text = " ".join([word.lower() for word in text.split(" ") if len(word) > 1])

        pid_element = root.find(".//article-id[@pub-id-type='pmid']")
        if pid_element is None:
            pid = ""
        else:
            pid = pid_element.text

        pmc_element = root.find(".//article-id[@pub-id-type='pmc']")
        if pmc_element is None:
            pmc = ""
        else:
            pmc = pmc_element.text

        doi_element = root.find(".//article-id[@pub-id-type='doi']")
        if doi_element is None:
            doi = ""
        else:
            doi = doi_element.text

        title_element = root.find(".//title-group/article-title")
        title = ET.tostring(title_element, encoding="utf-8", method="text").decode(
            "utf-8"
        )

        journal = root.find(".//journal-title-group/journal-title").text
        journal_id = root.find(".//journal-id[@journal-id-type='nlm-journal-id']").text

        issns = []

        journal_issnp_element = root.find(".//issn[@pub-type='ppub']")
        if journal_issnp_element is None:
            pass
        else:
            issns.append(journal_issnp_element.text)

        journal_issne_element = root.find(".//issn[@pub-type='epub']")
        if journal_issne_element is None:
            pass
        else:
            issns.append(journal_issne_element.text)

        year_element = root.find(".//pub-date[@pub-type='epub']/year")
        if year_element is None:
            year_element = root.find(".//pub-date[@pub-type='ppub']/year")

        if year_element is None:
            year = -1
        else:
            year = year_element.text

        doi = re.sub('"', "", doi)
        title = re.sub('"', "", title)
        # For some weird reason, some new-lines are not being replaced ...
        title = re.sub("\s\s+", " ", title).replace("\n", "")
        journal = re.sub('"', "", journal)
        issn = ",".join(issns)

        return (
            pid,
            pmc,
            '"' + doi + '"',
            year,
            '"' + title + '"',
            '"' + journal + '"',
            '"' + issn + '"',
            journal_id,
            '"' + text + '"',
        )


def main():
    with open("data.csv", "w+") as f_out:
        i = 0
        for path, _, files in os.walk("data"):
            for name in files:
                if i % 1000 == 0:
                    print(i)

                file_name = pathlib.PurePath(path, name)
                try:
                    f_out.write(",".join(process_article(file_name)) + "\n")
                except Exception as err:
                    print(err)
                i += 1


if __name__ == "__main__":
    main()
