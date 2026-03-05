import pandas as pd

# Leer archivo TSV
df = pd.read_excel("data/counts.xlsx")

df = df.iloc[:, 0:]
# Guardar como CSV
df.to_csv("data/counts.csv", index=False)