#%%
# Import a packages and make sure you have your compass-fiat environment activated
from fiat.io import open_csv
from pathlib import Path
import pandas as pd

#%%
# Load the CSV file into a pandas DataFrame
csv_path = Path(r"../computations/example/output", "output.csv")
df = pd.read_csv(csv_path)

# Check the column names to confirm the structure
print(df.columns)

# Extract the total damage values
total_damage = df['total_damage'].sum()
print(total_damage)

# %%
