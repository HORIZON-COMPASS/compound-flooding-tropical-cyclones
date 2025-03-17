#%%
# Sample Data
categories = ['Flood drivers']
data = {
    'Precipitation': [np.random.normal(loc=5, scale=2, size=100)],
    'Wind': [np.random.normal(loc=2, scale=1, size=100)],
    'SLR': [np.random.normal(loc=3, scale=1, size=100)],
    'Compound': [np.random.normal(loc=6, scale=1, size=100)],
}

# Define spacing parameters
n_sets = len(data)
n_categories = len(categories)
group_width = 0.8  # Reduce group width to tighten spacing
set_spacing = 0.2  # Less spacing between boxes
box_width = 0.12  # Keep boxes narrow
positions = []

# Compute positions for each box
for i, category in enumerate(categories):
    base = i * (group_width + 1.5)  # Adjust category spacing
    positions.extend([base + j * set_spacing for j in range(n_sets)])

# Flatten data for plotting
flattened_data = [item for sublist in zip(*data.values()) for item in sublist]

# Plot
fig, ax = plt.subplots(figsize=(8, 5))  # Smaller figure size
box = ax.boxplot(flattened_data, positions=positions, widths=box_width, patch_artist=True)

# Add colors for each set
colors = ['lightblue', 'lightgreen', 'coral', 'pink']
for patch, color in zip(box['boxes'], colors * n_categories):
    patch.set_facecolor(color)

# Custom x-ticks for categories
category_centers = [(i * (group_width + 1.5) + (group_width - set_spacing) / 2) for i in range(n_categories)]
ax.set_xticks(category_centers)
ax.set_xticklabels(categories, fontsize=14)

# Labels and grid
ax.set_ylabel('Diff in flood volume [%]', fontsize=14)
ax.set_title('Climate Attribution of Driver Contributions to Compound Flood from TC Idai', fontsize=16)
ax.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Ensure x-axis label is properly centered
ax.set_xlabel('')

# Legend
handles = [plt.Line2D([0], [0], color=color, lw=4, label=key) for key, color in zip(data.keys(), colors)]
ax.legend(handles, data.keys(), title='Driver', loc='upper left', fontsize=14, title_fontsize=14)

plt.tight_layout()
plt.show()


#%%
# Define flood types and their suffixes
flood_types = {
    'Pluvial': 'pluv',
    'Fluvial': 'fluv',
    'Coastal': 'coast',
    'Compound': 'cmpd'
}

# Define drivers
drivers = ['precipitation', 'temperature', 'wind', 'noSLR']

# Create dictionaries to store values for each flood type and driver
driver_contributions = {
    'precipitation': {},
    'temperature': {},
    'wind': {},
    'noSLR': {}
}

# Iterate through flood types and drivers to generate np.array datasets
for flood_type, suffix in flood_types.items():
    for driver in drivers:
        # Filter and create an array for the current flood type and driver
        driver_values = np.array([
            d['Volume_diff_from_F(%)'] for d in flood_characs
            if d['flood_type'] == flood_type and d['scenario'] == 'Counterfactual' and d['driver'] == driver
        ])
        
        # Ensure an empty array is created if no matching values are found
        if driver_values.size == 0:
            driver_values = np.array([])

        # Use the flood type suffix as the key
        dataset_name = f"{suffix}"
        driver_contributions[driver][dataset_name] = driver_values

#%%
### Make the plot ###
categories = ['Pluvial', 'Fluvial', 'Coastal', 'Compound']
data = {
    'Precipitation': [
        driver_contributions['precipitation']['pluv'],
        driver_contributions['precipitation']['fluv'],
        driver_contributions['precipitation']['coast'],
        driver_contributions['precipitation']['cmpd']
    ],
    'Temperature': [
        driver_contributions['temperature']['pluv'],
        driver_contributions['temperature']['fluv'],
        driver_contributions['temperature']['coast'],
        driver_contributions['temperature']['cmpd']
    ],
    'Wind': [
        driver_contributions['wind']['pluv'],
        driver_contributions['wind']['fluv'],
        driver_contributions['wind']['coast'],
        driver_contributions['wind']['cmpd']
    ],
    'SLR': [
        driver_contributions['noSLR']['pluv'],
        driver_contributions['noSLR']['fluv'],
        driver_contributions['noSLR']['coast'],
        driver_contributions['noSLR']['cmpd']
    ],
}

# Prepare the data for seaborn
flat_data = []
for category, values in zip(categories, zip(*data.values())):
    for driver, driver_values in zip(data.keys(), values):
        flat_data.extend([(category, driver, value) for value in driver_values])

df = pd.DataFrame(flat_data, columns=['Category', 'Driver', 'Value'])

# Plot
plt.figure(figsize=(12, 6))
sns.swarmplot(data=df, x='Category', y='Value', hue='Driver', dodge=True, palette='Set2')

# Customization
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
plt.title('Climate Attribution of Driver Contributions to Compound Flood from TC Idai')
plt.ylabel('Diff in Flood Volume [%]')
plt.xlabel('Flood Type')
plt.legend(title='Driver Contributions', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()





# %% 
                      #####################
                    ###### EXAMPLES ########
                      #####################
# Sample data
categories = ['Pluvial', 'Fluvial', 'Coastal', 'Compound']
data = {
    'Pluvial': [np.random.normal(loc=5, scale=2, size=5),
                np.random.normal(loc=4, scale=2, size=5),
                np.random.normal(loc=6, scale=2, size=5),
                np.random.normal(loc=6, scale=2, size=5)],
    'Fluvial': [np.random.normal(loc=-3, scale=1.5, size=5),
                np.random.normal(loc=-4, scale=1.5, size=5),
                np.random.normal(loc=-2, scale=1.5, size=5),
                np.random.normal(loc=-2, scale=1.5, size=5)],
    'Coastal': [np.random.normal(loc=2, scale=1, size=5),
                np.random.normal(loc=3, scale=1, size=5),
                np.random.normal(loc=1, scale=1, size=5),
                np.random.normal(loc=1, scale=1, size=5)],
    'Compound': [np.random.normal(loc=2, scale=1, size=5),
                 np.random.normal(loc=3, scale=1, size=5),
                 np.random.normal(loc=1, scale=1, size=5),
                 np.random.normal(loc=1, scale=1, size=5)],
}

# Create positions for the boxplots
n_sets = len(data)
n_categories = len(categories)
group_width = 1  # Width of each group of boxes
set_spacing = group_width / n_sets
positions = []

for i, category in enumerate(categories):
    base = i * (group_width + 0.5)  # 0.5 adds spacing between categories
    positions.extend([base + j * set_spacing for j in range(n_sets)])

# Flatten the data for plotting
flattened_data = [item for sublist in zip(*data.values()) for item in sublist]

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
box = ax.boxplot(flattened_data, positions=positions, widths=0.2, patch_artist=True)

# Add colors for each set
colors = ['lightblue', 'lightgreen', 'coral', 'gold']  # Extended to include 4 colors
for patch, color in zip(box['boxes'], colors * n_categories):
    patch.set_facecolor(color)

# Custom x-ticks for categories
category_centers = [(i * (group_width + 0.5) + (group_width - set_spacing) / 2) for i in range(n_categories)]
ax.set_xticks(category_centers)
ax.set_xticklabels(categories)

# Labels and grid
ax.set_ylabel('Diff in Flood Volume [%]')
ax.set_title('Boxplot with Grouped Boxes and Spaced Categories')
ax.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
ax.grid(axis='y', linestyle='--', alpha=0.7)

# Legend
handles = [plt.Line2D([0], [0], color=color, lw=4, label=key) for key, color in zip(data.keys(), colors)]
ax.legend(handles, data.keys(), title='Data Sets', loc='upper left')

plt.tight_layout()
plt.show()


# %%
# MAKE A VIOLIN PLOT
# Prepare the data for seaborn
flat_data = []
for category, values in zip(categories, zip(*data.values())):
    for dataset_index, dataset_values in enumerate(values):
        flat_data.extend([(category, f"Dataset {dataset_index+1}", value) for value in dataset_values])

df = pd.DataFrame(flat_data, columns=['Category', 'Dataset', 'Value'])

# Plot
plt.figure(figsize=(12, 6))
sns.violinplot(data=df, x='Category', y='Value', hue='Dataset', dodge=True, palette='Set2', inner='point')

# Customization
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
plt.title('Climate Attribution Data by Category (Violin Plot)')
plt.ylabel('Diff in Flood Volume [%]')
plt.xlabel('Flood Type')
plt.legend(title='Driver', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()

#%%

# Sample data for drivers (precipitation, temperature, etc.)
categories = ['Pluvial', 'Fluvial', 'Coastal', 'Compound']
data = {
    'Pluvial': [np.random.normal(loc=5, scale=2, size=10),
                np.random.normal(loc=4, scale=2, size=10),
                np.random.normal(loc=6, scale=2, size=10),
                np.random.normal(loc=6, scale=2, size=10)],
    'Fluvial': [np.random.normal(loc=-3, scale=1.5, size=10),
                np.random.normal(loc=-4, scale=1.5, size=10),
                np.random.normal(loc=-2, scale=1.5, size=10),
                np.random.normal(loc=-2, scale=1.5, size=10)],
    'Coastal': [np.random.normal(loc=2, scale=1, size=10),
                np.random.normal(loc=3, scale=1, size=10),
                np.random.normal(loc=1, scale=1, size=10),
                np.random.normal(loc=1, scale=1, size=10)],
    'Compound': [np.random.normal(loc=2, scale=1, size=10),
                 np.random.normal(loc=3, scale=1, size=10),
                 np.random.normal(loc=1, scale=1, size=10),
                 np.random.normal(loc=1, scale=1, size=10)],
}

# Extra base scenario data for the Pluvial Precipitation
extra_precipitation_data = np.random.normal(loc=7, scale=2, size=100)  # Extra data for Pluvial + Precipitation

# Prepare the data
flat_data = []
for category, values in zip(categories, zip(*data.values())):
    for dataset_index, dataset_values in enumerate(values):
        if category == 'Pluvial' and dataset_index == 0:  # Targeting Pluvial + Precipitation
            # Add the original and the extra data for Pluvial and Precipitation
            flat_data.extend([(category, f"Precipitation {dataset_index + 1}", value) for value in dataset_values])
            flat_data.extend([(category, f"Precipitation Extra", value) for value in extra_precipitation_data])
        else:
            flat_data.extend([(category, f"Dataset {dataset_index + 1}", value) for value in dataset_values])

df = pd.DataFrame(flat_data, columns=['Category', 'Dataset', 'Value'])

# Plot
plt.figure(figsize=(12, 6))

# Plot the regular violin plots for all categories and drivers
sns.violinplot(data=df, x='Category', y='Value', hue='Dataset', dodge=True, palette='Set2', inner='point', scale='count')

# Overlay the extra violin plot only for Pluvial Precipitation
extra_data = df[df['Dataset'] == 'Precipitation Extra']
sns.violinplot(data=extra_data, x='Category', y='Value', dodge=True, color='lightgrey', inner='point', scale='count', alpha=0.5)

# Customization
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Zero line for positive/negative values
plt.title('Climate Attribution Data by Category (Violin Plot) with Extra Precipitation Data')
plt.ylabel('Diff in Flood Volume [%]')
plt.xlabel('Flood Type')

# Adjust plot appearance
plt.legend(title='Driver', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()


#%%
# Example data
data = [2.5, 3.0, 2.8, 3.2, 2.9]

sns.swarmplot(data=[data])
plt.ylabel("Diff in Flood Volume [%]")
plt.title("Bee Swarm Plot for Small Sample Size")
plt.show()

# %%
sns.violinplot(data=[data], inner=None, color="lightblue")
sns.swarmplot(data=[data], color="black", size=7)
plt.ylabel("Values")
plt.title("Violin Plot with Data Points")
plt.show()
# %%
