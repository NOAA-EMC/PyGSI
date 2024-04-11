from SpatialTemporalStats import SpatialTemporalStats

# Set input and output paths
input_path = "/PATH/TO/Input/Files"
output_path = r'./Results'


# Set sensor name
sensor = "atms_n20"

# Set variable name and channel number
var_name = "Obs_Minus_Forecast_adjusted"
channel_no = 1

# Set start and end dates
start_date, end_date = '2023-03-01', '2023-03-10'

# Set region
# 1: global, 2: polar region, 3: mid-latitudes region,
# 4: tropics region, 5:southern mid-latitudes region, 6: southern polar region
region = 3

# Initialize SpatialTemporalStats object
my_tool = SpatialTemporalStats()

# Set resolution for grid generation
resolution = 2

# Generate grid
my_tool.generate_grid(resolution)  # Call generate_grid method)

# Set QC filter
QC_filter = True  # should be always False or true

# Set filter by variables
# can be an empty list
filter_by_vars=[]

#filter_by_vars = [("Land_Fraction", "lt", 0.9),]
# list each case in a separate tuple inside this list.
# options are 'lt' or 'gt' for 'less than' and 'greater than'

# Read observational values and perform analysis
o_minus_f_gdf = my_tool.read_obs_values(
    input_path,
    sensor,
    var_name,
    channel_no,
    start_date,
    end_date,
    filter_by_vars,
    QC_filter,
)

# Can save the results in a gpkg file
# o_minus_f_gdf.to_file("filename.gpkg", driver='GPKG')

# Plot observations
my_tool.plot_obs(o_minus_f_gdf, var_name, region, resolution, output_path)

# Make summary plots
summary_results = my_tool.make_summary_plots(
    input_path, sensor, var_name, start_date, end_date, QC_filter, output_path
)

# Print summary results
print(summary_results)
