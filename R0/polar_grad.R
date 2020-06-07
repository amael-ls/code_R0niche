
#### Aim of prog: plot the averaged gradients (north and south zones) on a radial plot

library(data.table)
library(tikzDevice)
library(stringi)
library(xtable)

options(max.print = 500)
rm(list = ls())

#### Tool functions
## Compute the azimuth of the averaged gradient. I assume dx != 0
# https://www.asprs.org/wp-content/uploads/pers/1987journal/aug/1987_aug_1109-1111.pdf
azimuth_grad = function(dx, dy, reverse = FALSE)
{
	if (reverse)
	{
		dx = -dx
		dy = -dy
	}
	return (ifelse(dx > 0, 90, 270) - 57.296*atan(dy/dx))
}

## Extract species from species_id
get_sp = function(species_id)
	return (stri_sub(str = species_id, from = stri_locate_first(str = species_id, regex = "-")[,1] + 1))

## Radial plot, assumed to be sorted by increasing angles for the labels
azimuthPlot = function(radius, theta, degree, species_labels, alternateArrows = TRUE, ...)
{
	# Transform degrees to radians
	if (degree)
		theta = theta*pi/180

	n = length(theta)

	if (length(radius) != 1 & length(radius) != n)
		print(paste0("Radius should be length one or length ", n, " (length of theta)"))

	# Azimuth => pi/2 is the origin, cos(pi/2 - a) = sin(a)
	x = radius*sin(theta)
	y = radius*cos(theta)

	if (alternateArrows)
	{
		if (length(radius == 1))
			radius = rep(radius, n)

		x[seq(1, n, 2)] = (radius[seq(1, n, 2)] + 0.3)*sin(theta[seq(1, n, 2)])
		y[seq(1, n, 2)] = (radius[seq(1, n, 2)] + 0.3)*cos(theta[seq(1, n, 2)])
	}

	# Add arrows
	sfsmisc::p.arrows(x1 = rep(0, n), y1 = rep(0, n), x2 = x, y2 = y, ...)

	# Add directions
	segments(0, -1, 0, 1, col = "#2A3344", lwd = 2)
	segments(-1, 0, 1, 0, col = "#2A3344", lwd = 2)
	text(x = 0, y = 1, labels = "N", pos = 3)
	text(x = 1, y = 0, labels = "E", pos = 4)
	text(x = 0, y = -1, labels = "S", pos = 1)
	text(x = -1, y = 0, labels = "W", pos = 2)

	# Add species' name
	pos_labels = integer(n)
	ind_west = x < 0
	text(x = x, y = y, labels = species_labels, pos = ifelse(ind_west, 2, 4), offset = 0.2)
}

#### Load data
## Gradient for s* = 0 and s* = 10
average_gradient = readRDS("average_gradient.rds")

## Compute the azimuth of the averaged gradient
average_gradient[, az_north_10 := azimuth_grad(dx = dx_north_10, dy = dy_north_10, reverse = TRUE)]
average_gradient[, az_south_10 := azimuth_grad(dx = dx_south_10, dy = dy_south_10, reverse = TRUE)]
average_gradient[, az_north_0 := azimuth_grad(dx = dx_north_0, dy = dy_north_0, reverse = TRUE)]
average_gradient[, az_south_0 := azimuth_grad(dx = dx_south_0, dy = dy_south_0, reverse = TRUE)]

## Extract species from species_id
average_gradient[, species := get_sp(species_id)]

#### Plots, canopy height s* = 10 m
# North
setorderv(average_gradient, "az_north_10")
tikz('./azimuth_north_10.tex', width = 6, height = 6)
op = par(mar = c(0, 0, 0, 0))
plot(0, pch = "", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), axes = FALSE,
	xlab = "", ylab = "", bg = "transparent")
azimuthPlot(1, average_gradient[, az_north_10], TRUE, average_gradient[, species], alternateArrows = TRUE, # seq(0.3, 1.6, 0.1)
	col = "#FD9859", fill = "#2058DC")
dev.off()

# South
setorderv(average_gradient, "az_south_10")
tikz('./azimuth_south_10.tex', width = 6, height = 6)
op = par(mar = c(0, 0, 0, 0))
plot(0, pch = "", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), axes = FALSE,
	xlab = "", ylab = "", bg = "transparent")
azimuthPlot(1, average_gradient[, az_south_10], TRUE, average_gradient[, species],
	col = "#FD9859", fill = "#2058DC")
dev.off()

## Save tex files
setorderv(average_gradient, "species")
gradient_xt = xtable(average_gradient[, .(species, dir_north_10, az_north_10, dir_south_10, az_south_10,
	dir_north_0, az_north_0, dir_south_0, az_south_0)], digits = 0)
print(gradient_xt, file = "./gradient.tex", include.rownames = FALSE, booktabs = TRUE)

## Plot all the species on different columns, s* = 10m
nbSpecies = average_gradient[, .N]
for (i in 1:nbSpecies)
{
	species = average_gradient[i, species]
	tikz(paste0("azimuth_", species, "_10.tex"), width = 1.5, height = 1.5)
	op = par(mar = c(0, 0, 0, 0))

	plot(x = NULL, y = NULL, xlim = c(-1.25, 1.25),
		ylim = c(-1.25, 1.25), axes = FALSE, xlab = "",
		ylab = "")

	## Centroids
	points(0, 0, pch = 15)
	points(0, 0.5, pch = 20)
	points(0, -0.5, pch = 20)

	## Compute and plot arrows
	# North
	theta = average_gradient[i, az_north_10]
	theta = theta*pi/180
	x = sin(theta)
	y = cos(theta)
	sfsmisc::p.arrows(x1 = 0, y1 = 0.5, x2 = x, y2 = 0.5 + y, col = "#2058DC", fill = "#2058DC")

	# South
	theta = average_gradient[i, az_south_10]
	theta = theta*pi/180
	x = sin(theta)
	y = cos(theta)
	sfsmisc::p.arrows(x1 = 0, y1 = -0.5, x2 = x, y2 = -0.5 + y, col = "#FD9859", fill = "#FD9859")

	## Add points (origin of arrows)
	points(0, 0.5, pch = 20, col = "#2058DC")
	points(0, -0.5, pch = 20, col = "#FD9859")

	dev.off()
}

#### Plots, canopy height s* = 0 m
# North
setorderv(average_gradient, "az_north_0")
tikz('./azimuth_north_0.tex', width = 6, height = 6)
op = par(mar = c(0, 0, 0, 0))
plot(0, pch = "", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), axes = FALSE,
	xlab = "", ylab = "", bg = "transparent")
azimuthPlot(1, average_gradient[, az_north_0], TRUE, average_gradient[, species], alternateArrows = TRUE, # seq(0.3, 1.6, 0.1)
	col = "#FD9859", fill = "#2058DC")
dev.off()

# South
setorderv(average_gradient, "az_south_0")
tikz('./azimuth_south_0.tex', width = 6, height = 6)
op = par(mar = c(0, 0, 0, 0))
plot(0, pch = "", xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), axes = FALSE,
	xlab = "", ylab = "", bg = "transparent")
azimuthPlot(1, average_gradient[, az_south_0], TRUE, average_gradient[, species],
	col = "#FD9859", fill = "#2058DC")
dev.off()

#### Plot all the species on different columns, s* = 0m
nbSpecies = average_gradient[, .N]
for (i in 1:nbSpecies)
{
	species = average_gradient[i, species]
	tikz(paste0("azimuth_", species, "_0.tex"), width = 1.5, height = 1.5)
	op = par(mar = c(0, 0, 0, 0))

	plot(x = NULL, y = NULL, xlim = c(-1.25, 1.25),
		ylim = c(-1.25, 1.25), axes = FALSE, xlab = "",
		ylab = "")

	## Centroids
	points(0, 0, pch = 15)
	points(0, 0.5, pch = 20)
	points(0, -0.5, pch = 20)

	## Compute and plot arrows
	# North
	theta = average_gradient[i, az_north_0]
	theta = theta*pi/180
	x = sin(theta)
	y = cos(theta)
	sfsmisc::p.arrows(x1 = 0, y1 = 0.5, x2 = x, y2 = 0.5 + y, col = "#2058DC", fill = "#2058DC")

	# South
	theta = average_gradient[i, az_south_0]
	theta = theta*pi/180
	x = sin(theta)
	y = cos(theta)
	sfsmisc::p.arrows(x1 = 0, y1 = -0.5, x2 = x, y2 = -0.5 + y, col = "#FD9859", fill = "#FD9859")

	## Add points (origin of arrows)
	points(0, 0.5, pch = 20, col = "#2058DC")
	points(0, -0.5, pch = 20, col = "#FD9859")

	dev.off()
}

## Legend to read graphs
tikz("legend_cols.tex", width = 1.5, height = 1.5)
op = par(mar = c(0, 0, 0, 0))

plot(x = NULL, y = NULL, xlim = c(-1.25, 1.25),
	ylim = c(-1.25, 1.25), axes = FALSE, xlab = "",
	ylab = "")

# Segments
segments(x0 = 0, y0 = -0.75, x1 = 0, y1 = 0.75, lwd = 0.15, lty = "dashed")
segments(x0 = -0.75, y0 = 0, x1 = 0.75, y1 = 0, lwd = 0.15, lty = "dashed")

# Centroids
points(0, 0, pch = 15)
points(0, 0.5, pch = 19, col = "#2058DC")
points(0, -0.5, pch = 19, col = "#FD9859")

# Compass directions (see remark at the beginning of this code)
text(0, 0.8, "N", pos = 3, offset = 0.25)
text(0.8, 0, "E", pos = 4, offset = 0.25)
text(0, -0.8, "S", pos = 1, offset = 0.25)
text(-0.8, 0, "W", pos = 2, offset = 0.25)

dev.off()
