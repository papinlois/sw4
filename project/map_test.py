# Kate code for the plot
fig = plt.figure(figsize=(12, 8))  # Increase the figure size for clarity

# Initial plot with all events and stations
ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
plt.scatter(stas['lon'], stas['lat'], marker='v', label='Stations')
[plt.text(i, j, f'{name}', va='top', ha='left') for (i, j, name) in zip(stas['lon'], stas['lat'], stas['name'])]
if multipleEvents:
    plt.scatter(eventLocs['lon'], eventLocs['lat'], marker='o', edgecolors='darkorange', facecolors='orange',
                label='Event Locations')

# Scatter plot for the event to model
plt.scatter(event['lon'], event['lat'], marker='*', facecolors='red', label="Event to Model")
# plt.text(event['lon'], event['lat'], f"({event['lat']:.3f},{event['lon']:.3f})")
plt.scatter(sta['lon'], sta['lat'], marker='*', facecolors='green', label="Choosen Station")

ax.legend()
ax.coastlines()
gl = ax.gridlines(draw_labels=True)
Plotlims = ax.get_extent()

# These values are roughly the rectangle that encompasses all of our points plus a little buffer
lonMin_ini, lonMax_ini, latMin_ini, latMax_ini = Plotlims
x_ini = [Plotlims[0], Plotlims[1], Plotlims[1], Plotlims[0], Plotlims[0]]
y_ini = [Plotlims[2], Plotlims[2], Plotlims[3], Plotlims[3], Plotlims[2]]

# Calculate the minimum longitude and latitude values for the zoomed-in map
lonMin_zoom, lonMax_zoom = min(sta['lon']), max(sta['lon'])
latMin_zoom, latMax_zoom = min(sta['lat']), max(sta['lat'])

# Add a little buffer around the zoomed-in region
buffer = 0.1
lonMin_zoom -= buffer
lonMax_zoom += buffer
latMin_zoom -= buffer
latMax_zoom += buffer

# Create a zoomed-in map
ax_zoom = fig.add_axes([0.6, 0.6, 0.25, 0.25], projection=ccrs.PlateCarree())  # [left, bottom, width, height]
ax_zoom.set_extent([lonMin_zoom, lonMax_zoom, latMin_zoom, latMax_zoom])
ax_zoom.scatter(stas['lon'], stas['lat'], marker='v', label='Stations')
[ax_zoom.text(i, j, f'{name}', va='top', ha='left') for (i, j, name) in zip(stas['lon'], stas['lat'], stas['name'])]
if multipleEvents:
    ax_zoom.scatter(eventLocs['lon'], eventLocs['lat'], marker='o', edgecolors='darkorange', facecolors='orange',
                     label='Event Locations')
ax_zoom.scatter(event['lon'], event['lat'], marker='*', facecolors='red', label="Event to Model")
ax_zoom.scatter(sta['lon'], sta['lat'], marker='*', facecolors='green', label="Choosen Station")
ax_zoom.gridlines(draw_labels=True)

plt.show()