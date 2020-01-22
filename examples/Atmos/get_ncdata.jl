# This example demonstrates extraction / usage of NETCDF data from 
# CFMIP archives to prescribe initial conditions within the CliMA 
# framework. Future examples will consider extensions using forcing 
# functions.

# Requires NCDatasets.jl Julia package
# Requires Plots.jl for visualisation

using NCDatasets
using Plots

# Load single dataset in read-only mode
data = Dataset("./cfsites_forcing.2010071518.nc","r");

# Load specific site group via numeric ID in NetCDF file (requires generalisation)
siteid= data.group["site22"];

# Allow strings to be read as varnames
# TODO: Is this the most robust solution? 
function str2var(str::String, var::Any)
  str = Symbol(str);
  @eval(($str)=($var));
end


# Load all variables
for (varname,var) in siteid
  str2var(varname,var[:,1]);
end

# List Variables (Optional Function)
function list_vars(siteid)
  for (varname,var) in siteid
    @show(varname)
  end
end

# Uncomment following line to list all variables for the current site-ID
# list_vars(siteid)

# Plot Attributes (Simple GKS)
fontfamily = "Arial"
fontsize = 6
layout = (2,2);
ylims = (0, 6);
lw= 2;

height = height / 1e3

p1 = plot(pfull, height, 
          xlabel = "Pressure[Pa]",
          ylabel = "Altitude[km]",
          ylims = ylims,
          linewidth = lw,
          label = "", 
          xtickfont = font(fontsize,fontfamily));
p2 = plot(temp, height, 
          xlabel = "Temperature[K]",
          ylabel = "Altitude[km]",
          ylims = ylims,
          linewidth = lw,
          label = "", 
          xtickfont = font(fontsize,fontfamily));
p3 = plot(sphum, height, 
          xlabel = "Spec. Humidity[1]",
          ylabel = "Altitude[km]",
          ylims = ylims,
          linewidth = lw,
          label = "", 
          xtickfont = font(fontsize,fontfamily));
p4 = plot(ucomp, height, 
          xlabel = "Windspeed[m/s]",
          ylabel = "Altitude[km]",
          ylims = ylims,
          linewidth = lw,
          label = "Westerly", 
          xtickfont = font(fontsize,fontfamily));
p4 = plot!(vcomp, height, 
           xlabel = "Windspeed[m/s]",
           ylabel = "Altitude[km]",
           ylims = ylims,
           linewidth = lw,
           label = "Southerly", 
           xtickfont = font(fontsize,fontfamily));

plot(p1,p2,p3,p4, layout=layout)
