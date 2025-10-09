#%%
import pandas as pd
from deep_translator import GoogleTranslator

# Load Excel
df = pd.read_excel("data/07. SOFALA .xlsx")

#%%
# Example: translate a column named "text"
df_translated = df.apply(lambda x: GoogleTranslator(source="auto", target="en").translate(str(x)))

# Save back to Excel
df_translated.to_excel("data/Census_Sofala_translated.xlsx", index=False)
# %%
import pandas as pd
from deep_translator import GoogleTranslator

# Load all sheets into a dict {sheet_name: DataFrame}
file_in = "data/07. SOFALA .xlsx"
file_out = "data/Census_Sofala_translated.xlsx"

sheets = pd.read_excel(file_in, sheet_name=None)

translator = GoogleTranslator(source='auto', target='en')

# Function to translate a whole DataFrame
def translate_df(df):
    return df.applymap(lambda x: translator.translate(x) if isinstance(x, str) else x)

# Translate all sheets
translated_sheets = {name: translate_df(df) for name, df in sheets.items()}

# Save back to Excel (all sheets preserved)
with pd.ExcelWriter(file_out, engine="openpyxl") as writer:
    for name, df in translated_sheets.items():
        df.to_excel(writer, sheet_name=name, index=False)

print(f"Translated Excel saved to {file_out}")

# %%
