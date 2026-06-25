import pypdf
import re

pdf_path = r"C:\Users\Nathan\Downloads\qcad050.pdf"
reader = pypdf.PdfReader(pdf_path)

print(f"Total Pages: {len(reader.pages)}")

# Let's search for some patterns
for idx, page in enumerate(reader.pages):
    text = page.extract_text()
    page_num = idx + 1
    
    # Check if this page contains Table 2 or cost values
    if "Table 2" in text or "1,899" in text or "1899" in text or "3,569" in text or "3569" in text:
        print(f"--- MATCH on Page {page_num} ---")
        lines = text.split("\n")
        for line in lines:
            if any(term in line for term in ["Table 2", "1,899", "1899", "3,569", "3569", "1,925", "1925", "Cost", "QALY"]):
                # Print safely replacing non-ascii characters
                safe_line = line.encode('ascii', errors='replace').decode('ascii')
                print(f"  {safe_line}")

# Let's dump all text to a markdown file for thorough inspection
with open("scratch/paper_text.txt", "w", encoding="utf-8") as f:
    for idx, page in enumerate(reader.pages):
        f.write(f"\n\n--- PAGE {idx+1} ---\n\n")
        f.write(page.extract_text())
print("Saved all text to scratch/paper_text.txt")
