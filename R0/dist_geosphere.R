# Author: George Wang & Robert J. Hijmans
# August 2010
# version 1
# license GPL3

my_alongTrackDistance <- function(p1, p2, p3, r=6378137) {
	toRad <- pi / 180
	p1 <- .pointsToMatrix(p1)
	p2 <- .pointsToMatrix(p2)
	p3 <- .pointsToMatrix(p3)
	p <- cbind(p1[,1], p1[,2], p2[,1], p2[,2], p3[,1], p3[,2], as.vector(r))
	p1 <- p[,1:2,drop=FALSE]
	p2 <- p[,3:4,drop=FALSE]
	p3 <- p[,5:6,drop=FALSE]
	r = p[,7]

	tc <- bearing(p1, p2) * toRad
	tcp <- bearing(p1, p3) * toRad
    dp <- distCosine(p1, p3, r=1)
	xtr <- asin(sin(tcp-tc) * sin(dp))

# +1/-1 for ahead/behind [lat1,lon1]
	bearing <- sign(cos(tc - tcp))

	vec <- cos(dp) / cos(xtr)

	# Add a security here, for numbers not in [-1, 1] but very close
	if (any(abs(vec) > 1))
	{
		l_sup <- length(vec[vec > 1]) != 0
		l_inf <- length(vec[vec < -1]) != 0

		if (l_sup & all.equal(vec[vec > 1], rep(1, length(vec[vec > 1]))))
			vec[vec > 1] <- 1
		if (l_inf & all.equal(vec[vec < -1], rep(-1, length(vec[vec < -1]))))
			vec[vec < -1] <- -1
	}

	dist <- bearing * acos(vec) * r

	if (is.vector(dist)) { dist <- matrix(dist) }
	colnames(dist) <- 'distance'

	return(abs(dist))
}


.pointsToMatrix <- function(p, checkLonLat=TRUE, poly=FALSE) {
	if (inherits(p, 'SpatialPoints')) {
		test <- !is.projected(p)
		if (! isTRUE (test) ) {
			if (is.na(test)) {
				warning('Coordinate reference system of SpatialPoints object is not set. Assuming it is degrees (longitude/latitude)!')
			} else {
				stop('Points are projected. They should be in degrees (longitude/latitude)')
			}
			# or rather transform them ....?
		}
		p <- coordinates(p)
	} else if (is.data.frame(p)) {
		p <- as.matrix(p)
	} else

	if (is.vector(p)){
		if (length(p) != 2) {
			stop('Wrong length for a vector, should be 2')
		} else {
			p <- matrix(p, ncol=2)
		}
	} else if (is.matrix(p)) {
		if (ncol(p) != 2) {
			stop( 'A points matrix should have 2 columns')
		}
		cn <- colnames(p)
		if (length(cn) == 2) {
			if (toupper(cn[1]) == 'Y' | toupper(cn[2]) == 'X')  {
				warning('Suspect column names (x and y reversed?)')
			}
			if (toupper(substr(cn[1],1,3) == 'LAT' | toupper(substr(cn[2],1,3)) == 'LON'))  {
				warning('Suspect column names (longitude and latitude reversed?)')
			}
		}
	} else {
		stop('points should be vectors of length 2, matrices with 2 columns, or inheriting from a SpatialPoints* object')
	}

	if (! is.numeric(p) ) { p[] <- as.numeric(p) }

	if (checkLonLat & nrow(p) > 0) {
		if (length(stats::na.omit(p[,1])) > 0) {
			if (min(p[,1], na.rm=TRUE) < -360) { stop('longitude < -360') }
			if (max(p[,1], na.rm=TRUE) > 360) {  stop('longitude > 360')  }
			if (min(p[,1], na.rm=TRUE) < -180) { warning('longitude < -180') }
			if (max(p[,1], na.rm=TRUE) > 180) {  warning('longitude > 180')  }
		}
		if (length(stats::na.omit(p[,2])) > 0) {
			if (min(p[,2], na.rm=TRUE) < -90) {  stop('latitude < -90')  }
			if (max(p[,2], na.rm=TRUE) > 90) {  stop('latitude > 90')  }
		}
	}


	if (poly) {
		if (! isTRUE(all.equal(p[1,], p[nrow(p),]))) {
			p <- rbind(p, p[1,])
		}

		i <- p[-nrow(p),1] == p[-1,1] &  p[-nrow(p),2] == p[-1,2]
		i <- which(isTRUE(i))
		if (length(i) > 0) {
			p <- p[-i, ,drop=FALSE]
		}

		.isPolygon(p)
	}

	return(p)
}


.spDistPoint2Line <- function(p, line, distfun) {
	test <- !is.projected(line)
	if (! isTRUE (test) ) {
		if (is.na(test)) {
			warning('Coordinate reference system of SpatialPolygons object is not set. Assuming it is degrees (longitude/latitude)!')
		} else {
			stop('Points are projected. They should be in degrees (longitude/latitude)')
		}
		# or rather transform them ....?
	}

	x <- line@lines
	n <- length(x)
	res <- matrix(nrow=nrow(p), ncol=4)
	colnames(res) <- c("distance","lon","lat","ID")
	res[] <- Inf

	for (i in 1:n) {
		parts <- length(x[[i]]@Lines )
		for (j in 1:parts) {
			crd <- x[[i]]@Lines[[j]]@coords
			r <- cbind(my_dist2Line(p, crd, distfun), i)
			k <- r[,1] < res[,1]
			res[k, ] <- r[k, ]
		}
	}
	return(res)
}


my_dist2Line <- function(p, line, distfun=distGeo) {

	p <- .pointsToMatrix(p)

	if (inherits(line, 'SpatialPolygons')) {
		line <- methods::as(line, 'SpatialLines')
	}
	if (inherits(line, 'SpatialLines')) {
		return( .spDistPoint2Line(p, line, distfun) )
	}

	line <- .pointsToMatrix(line)
	line1 <- line[-nrow(line), ,drop=FALSE]
	line2 <- line[-1, ,drop=FALSE]
	seglength  <- distfun(line1, line2)

	res <- matrix(nrow=nrow(p), ncol=3)
	colnames(res) <- c("distance","lon","lat")

	for (i in 1:nrow(p)) {
		xy <- p[i,]
# the shortest distance of a point to a great circle
		crossdist <- abs(dist2gc(line1, line2, xy))

# the my_alongTrackDistance is the length of the path along the great circle to the point of intersection
# there are two, depending on which node you start
# we want to use the min, but the max needs to be < segment length
		trackdist1 <- my_alongTrackDistance(line1, line2, xy)
		trackdist2 <- my_alongTrackDistance(line2, line1, xy)
		mintrackdist <- pmin(trackdist1, trackdist2)
		maxtrackdist <- pmax(trackdist1, trackdist2)
		crossdist[maxtrackdist >= seglength] <- NA

# if the crossdist is NA, we use the distance to the nodes
		nodedist <- distfun(xy, line)

		warnopt = getOption('warn')
	 	options('warn'=-1)
		distmin1 <- min(nodedist, na.rm=TRUE)
		distmin2 <- min(crossdist, na.rm=TRUE)
		options('warn'= warnopt)

		if (distmin1 <= distmin2) {
			j <- which.min(nodedist)
			res[i,] <- c(distmin1, line[j,])
		} else {
			j <- which.min(crossdist)
			# if else to determine from which node to start
			if (trackdist1[j] < trackdist2[j]) {
				bear <- bearing(line1[j,], line2[j,])
				pt <- destPoint(line1[j,], bear, mintrackdist[j])
				res[i,] <- c(crossdist[j], pt)
			} else {
				bear <- bearing(line2[j,], line1[j,])
				pt <- destPoint(line2[j,], bear, mintrackdist[j])
				res[i,] <- c(crossdist[j], pt)
			}
		}
	}
	return(res)
}
