library(rgdal)
library(plotrix)
library(rgeos)
library(data.table)
library(magrittr)

sites <- fread('../data/combined/all_samples.samples.meta.csv') %>%
         .[, .(Lat = unique(latitude), Long = unique(longitude)), by = c('location', 'country')]
sites[, species := ifelse(location %in% c('Avrankou', 'Korle-Bu'), 'coluzzii', 'gambiae')]

pal <- c(gambiae = rgb(0.101, 0.301, 0.506), 
         coluzzii = rgb(0.89, 0.102, 0.11))

# Load the map of Africa
Afmap <- readOGR(dsn = 'Africa_SHP', layer = 'Africa')

sites[, ':='(label.offset.lat = 0.45,
             label.offset.long = -0.1)
]
sites[location %in% c('Avrankou', 'Korle-Bu'), ':='(label.offset.lat = -0.45, label.offset.long = 0.2)]

shadow.text <- function(x, y, label, shadow.size, shadow.col = 'white', text.col = 'black', ...){
	for (i in c(1,1,-1,-1)){
		for (j in c(1,-1,1,-1)){
			 text(x + shadow.size*i, y + shadow.size*j, label, col = shadow.col, ...)
		}
	}
	text(x, y, label, col = text.col, ...)
}

lighten.col <- function(color, lightness, alpha = alpha){
	col.rgb <- col2rgb(color)/255
	rgb(t(1-(1-col.rgb)*lightness), alpha = alpha)
}

# Plot the map 
plot.collection.sites <- function(sites, palette){
	sites <- copy(sites)
	par(mar = c(0.5, 0.5, 0.5, 0.5))
	countries <-  unique(sites$country) %>%
	              {setNames(sub("'", "`", .), .)}
	plot(Afmap[Afmap$COUNTRY %in% countries,], border = 'black', col = 'white', lwd = 2, bg = rgb(0.482, 0.694, 0.741))
	# Make the other countries grey
	plot(Afmap[!(Afmap$COUNTRY %in% countries),], border = 'grey60', col = 'grey80', add = T)
	# Draw the main countries again to the borders are bold all the way
	plot(Afmap[Afmap$COUNTRY %in% countries,], border = 'black', col = 'white', lwd = 2, add = T)
	# Outline the focal countries in bold
	# Label the countires
	for (country in names(countries)){
		centroid = gCentroid(Afmap[Afmap$COUNTRY == countries[country], ])
		shadow.text(centroid$x, centroid$y+0.5, country, 0.02, text.col = 'grey20', cex = 1.3, font = 2)
	}
	# Add the GAARDIAN sampling locations. 
	# Tiassale
	sites[, 
		{draw.circle(Long, Lat, radius = 0.1, col = lighten.col(palette[species], 0.5), nv = 1000)
	 	 shadow.text(Long+label.offset.long, Lat+label.offset.lat, text.col = palette[species], location, 0.02, cex = 1.2)},
		by = 'location'
	]
}

#pdf('collection_sites.pdf', width = 7, height = 4)
#plot.collection.sites(sites, pal)
#dev.off()

png('collection_sites.png', width = 7, height = 4, units = 'in', res = 300)
plot.collection.sites(sites, pal)
dev.off()


