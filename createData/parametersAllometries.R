
#### Aim of prog:
## Calculate the height of a tree given its dbh
# 		- Parameters for 106 species, source table S3, Purves2007, a and b
# 		- Function define at the end of the program, source appendix S2, end 1st page
## Calculate the radius of the crown at a distance y from the top, given dbh of a tree
#		- Parameters for 106 species, source table S3, Purves2007, T
#		- Function define ate the end of the program, source appendix S1, eq S1.6, S1.7
## Define the C0-C1 data table
#		- Source = table S2 Purves2007
#
# The two following species have synonyms in the database:
#		TAX-ASC: Taxodium ascendens
#		https://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=183433#null
#
#		NA-CAR-ALB: Carya alba
#		https://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=501306#null
#
## Bibliography:
#		- Crown plasticity and competition for canopy space: a new spatially implicit model parameterized for 250 North American tree species
#
#		- Managing understory light conditions in boreal mixedwoods through variation in the intensity and spatial pattern of harvest: A modelling approach
#
#		- Predictions of understorey light conditions in northern hardwood forests following parameterization, sensitivity analysis, and tests of the SORTIE light model

library(data.table)
library(stringi)

#### Tool function
cleanName = function(str)
	return(stri_sub(str, from = stri_locate_first(str, regex = "-")[1] + 1))

#### Create the data table containing the parameters from Purves 2007
## Species names
ls_species = c("18032-ABI-BAL", "18034-PIC-RUB", "28728-ACE-RUB", "32931-FRA-AME",
"28754-ACE-PEN", "27806-COR-FLO", "19481-BET-ALL", "19337-QUE-GAM", "19408-QUE-RUB",
"183400-TSU-HET", "183327-PIN-CON", "195773-POP-TRE", "19280-QUE-NIG", "183397-TSU-CAN", "19020-PLA-OCC",
"28731-ACE-SAC", "18070-MAG-VIR", "NA-TAX-ASC", "19487-BET-LEN", "19027-LIQ-STY", "501145-CAL-DEC", "19447-QUE-VEL",
"183335-PIN-ECH", "183302-PIC-MAR", "18158-SAS-ALB", "18086-LIR-TUL", "32929-FRA-PEN", "19366-QUE-KEL", "NA-CAR-ALB",
"23855-DIO-VIR", "19050-ULM-RUB", "18048-JUN-VIR", "19246-CAR-TEX", "183394-PIN-VIR", "19051-ULM-ALA", "19464-LIT-DEN",
"183424-PSE-MEN", "19511-OST-VIR", "27982-ILE-OPA", "19040-CEL-OCC", "19462-FAG-GRA", "19287-QUE-MAC",
"183336-PIN-EDU", "19327-QUE-ELL", "19049-ULM-AME", "19283-QUE-VIR", "28757-ACE-SAC", "18037-PIN-TAE",
"183284-ABI-GRA", "27821-NYS-SYL", "32945-FRA-NIG", "23690-OXY-ARB", "19290-QUE-ALB", "183385-PIN-STR",
"18044-THU-PLI", "19277-QUE-FAL" , "19489-BET-PAP", "19368-QUE-LAU", "504804-ROB-PSE", "19227-CAR-COR",
"19504-CAR-CAR", "NA-QUE-PRI", "19254-JUG-NIG", "19042-CEL-LAE", "183417-LAR-OCC", "183319-PIN-BAN", "18038-PIN-PAL",
"183291-PIC-ENG", "505490-THU-OCC", "22463-POP-GRA", "19474-ALN-RUB", "22484-SAL-NIG", "19242-CAR-OVA",
"183412-LAR-LAR", "19422-QUE-STE", "28749-ACE-NEG" , "183365-PIN-PON", "19312-QUE-CHR", "195925-NYS-BIF",
"183295-PIC-GLA", "19231-CAR-GLA", "19288-QUE-COC", "18041-TAX-DIS", "183311-PIN-ALB", "181826-ABI-CON",
"19497-BET-POP", "NA-CHA-NOO", "18036-PIN-ELL", "183309-PIC-SIT", "19374-QUE-MAR", "24799-PRU-PEN", "183375-PIN-RES",
"19282-QUE-PHE", "194855-JUN-OCC", "26879-PRO-GLA", "22453-POP-BAL", "183353-PIN-MON", "194812-JUN-ASH",
"25134-CER-LED", "194872-JUN-SCO", "183345-PIN-JEF", "194859-JUN-OST", "27824-NYS-AQU", "181834-ABI-MAG",
"21536-TIL-AME", "24764-PRU-SER", "181824-ABI-AMA", "181830-ABI-LAS", "183402-TSU-MER")

## Intercept a dbh--height allometry
ls_a = c(.1203, .1309, .4769, .4466,
	.4352, .4421, .5297, .3632, .4363,
	.1829, .1645, .351, .4053, .2893, .4994,
	.4742, .4052, .198, .492, .3699, .05112, .3725,
	.2356, .175, .4113, .4368, .4525, .3228, .3633,
	.3575, .394, .3195, .3222, .2829, .3401, .3184,
	.1426, .4454, .2903, .3411, .3989, .3424,
	.1773, .4155, .4192, .3613, .4801, .2401,
	.04658, .3566, .4782, .3931, .3989, .228,
	.2029, .3558, .4423, .4013, .4268, .3579,
	.5218, .4206, .3885, .423, .316, .223, .2001,
	.06969, .2302, .4359, .392, .4914, .3742,
	.3045, .3291, .5276, .06788, .2034, .464,
	.1404, .3502, .4052, .1858, .09986, .03303,
	.4551, .1417, .1785, .1663, .3371, .4275, .2567,
	.3963, .2167, .3856, .4062, .04624, .2366,
	.4058, .2368, .1773, .2369, .4224, .0511,
	.4258, .4527, -0.09029, .01819, -0.05736)

## Coefficient b dbh--height allometry
ls_b = c(.783, .7555, .5641, .5645,
	.5372, .4583, .471, .7387, .5678,
	.7327, .77, .6078, .6223, .6058, .5417,
	.5743, .6132, .7415, .5399, .6703, .7242, .6039,
	.7331, .7487, .5707, .6302, .534, .5092, .6506,
	.6226, .6015, .5537, .6448, .691, .6205, .6164,
	.7589, .5438, .6099, .6015, .5954, .5546,
	.6372, .5438, .5523, .5469, .5341, .7302,
	.7875, .6207, .524, .5876, .6029, .6755,
	.6842, .631, .5634, .5799, .5418, .6582,
	.4359, .5805, .5843, .5687, .7169, .6828, .7659,
	.7796, .6194, .6066, .6299, .4416, .6423,
	.6559, .6147, .4269, .7441, .5863, .5688,
	.7354, .6749, .5969, .7323, .708, .7735,
	.5319, .6981, .7925, .724, .5355, .5369, .6522,
	.5934, .5037, .5051, .5705, .8043, .5478,
	.5224, .5566, .6977, .4802, .5739, .7997,
	.5768, .525, .8955, .8281, .8255)

## Parameters T for height--crown allometry
ls_T = c(.278, .324, .536, .551,
	.674, .624, .592, .359, .538,
	.337, .235, .404, .590, .464, .578,
	.560, .454, .623, .571, .432, .240, .507,
	.407, .201, .450, .489, .511, .457, .599,
	.526, .603, .380, .686, .427, .717, .365,
	.335, .557, .462, .599, .626, .489,
	NA, .508, .536, .639, .486, .355,
	.292, .503, .374, .475, .543, .417,
	.334, .639, .435, .634, .512, .577,
	.719, .509, .653, .670, .344, .342, .462,
	.208, .251, .424, .437, .523, .486,
	.308, .526, .646, .314, .349, .470,
	.278, .564, .576, .361, .290, .236,
	.371, .202, .406, .316, .380, .482, .287,
	.579, .015, .635, .345, .172, NA,
	.222, .716, .292, NA, .485, .181,
	.462, .520, .264, .181, .209)

#### Data table
## Creating
purves2007_allometries = data.table(species = ls_species, a = ls_a, b = ls_b, T = ls_T)

## Discard TSN from species code
purves2007_allometries[, sp := cleanName(species), by = species]

## Sort by species name
setorderv(purves2007_allometries, cols = "sp")

## Save the allometries to check with data table from readS3_purves2007.R
saveRDS(purves2007_allometries, "./purves2007_allometries.rds")

#### Define allometries
## From dbh to height, appendix S2
dbhToHeight = function(dbh, a, b, mm = TRUE) # If in mm, then dbh/10, i.e., 10^{-1}
	return (10^(a - mm*b + b*log10(dbh)))

## From heigth to dbh, appendix S2 /!\ Return dbh in cm /!\
heightToDbh = function(height, a, b)
	return (ifelse(height > 0, 10^((log10(height) - a) / b), 0))

## From height to crown radius at a given distance to the top (using dbh as an intermediate), appendix S1 and S3. Height in m, distToTop is the distance to the top at which the crown radius is evaluated. Finally a, b, T_param and C0_C1 are allometries from Purves 2007. I need to use M[i] and B[i] here, as I give vectors fo heights of different species (M and B are sp-specific due to T_param)
heightToCrownArea = function(height, distToTop, a, b, T_param, C0_C1)
{
	# Table S2
	R0_C0 = C0_C1[parameter == "R0", C0]
	R0_C1 = C0_C1[parameter == "R0", C1]

	R40_C0 = C0_C1[parameter == "R40", C0]
	R40_C1 = C0_C1[parameter == "R40", C1]

	M_C0 = C0_C1[parameter == "M", C0]
	M_C1 = C0_C1[parameter == "M", C1]

	B_C0 = C0_C1[parameter == "B", C0]
	B_C1 = C0_C1[parameter == "B", C1]

	# Appendix S3, Eq S3.3 (erroneously denoted S2.3 in the article)
	R0 = (1 - T_param)*R0_C0 + T_param*R0_C1
	R40 = (1 - T_param)*R40_C0 + T_param*R40_C1

	M = (1 - T_param)*M_C0 + T_param*M_C1
	B = (1 - T_param)*B_C0 + T_param*B_C1

	# All the following part is vectorised
	# Convert height to dbh. /!\ If height is in m, then dbh is in cm /!\
	dbh = heightToDbh(height, a, b)

	# Calculate potential max radius Eq S1.6, if height = 0, then dbh = 0
	Rp_max = R0 + (R40 - R0)*dbh/40

	n = length(height)
	crownRadius = 0*numeric(length = n) # To be sure it is initialised to zero

	# Vector crown radius
	for (i in 1:n)
		if (height[i] > 0)
			crownRadius[i] = Rp_max[i] * ( min(distToTop[i], height[i]*M[i]) / (height[i]*M[i]) )^B[i] # R_{i, y}^p

	crownArea = pi * crownRadius * crownRadius

	return (crownArea)
}

## From height to crown radius at a given distance to the top (using dbh as an intermediate), appendix S1 and S3. Dbh is in mm, phi_star is the s* (for diameter). Finally a, b, T_param and C0_C1 are allometries from Purves 2007
dbhToCrownArea = function(dbh, phi_star, a, b, T_param, C0_C1, mm = TRUE)
{
	dbh[dbh < 0] = 0 # Due to precision, I get sometimes negative close to zero values
	phi_star[phi_star < 0] = 0 # Due to precision, I get sometimes negative close to zero values

	# Table S2
	R0_C0 = C0_C1[parameter == "R0", C0]
	R0_C1 = C0_C1[parameter == "R0", C1]

	R40_C0 = C0_C1[parameter == "R40", C0]
	R40_C1 = C0_C1[parameter == "R40", C1]

	M_C0 = C0_C1[parameter == "M", C0]
	M_C1 = C0_C1[parameter == "M", C1]

	B_C0 = C0_C1[parameter == "B", C0]
	B_C1 = C0_C1[parameter == "B", C1]

	# Appendix S3, Eq S3.3 (erroneously denoted S2.3 in the article)
	R0 = (1 - T_param)*R0_C0 + T_param*R0_C1
	R40 = (1 - T_param)*R40_C0 + T_param*R40_C1

	M = (1 - T_param)*M_C0 + T_param*M_C1
	B = (1 - T_param)*B_C0 + T_param*B_C1

	# Calculate potential max radius Eq S1.6, /!\ dbh in cm /!\
	Rp_max = R0 + (R40 - R0)*dbh/40
	if (mm)
		Rp_max = R0 + (R40 - R0)*dbh/400

	# All the following part is vectorised
	# Convert dbh to height. /!\ If dbh is in mm, then height is in m /!\
	height = dbhToHeight(dbh, a, b)
	height_star = dbhToHeight(phi_star, a, b)

	# Calculate distance to the top
	distToTop = height - height_star

	print(B)
	print(M)
	print(distToTop)

	n = length(height)
	crownRadius = 0*numeric(length = n) # To be sure it is initialised to zero
	# Vector crown radius
	for (i in 1:n)
		if (height[i] > 0)
			crownRadius[i] = Rp_max[i] *
			( min(distToTop[i], height[i]*M) / (height[i]*M) )^B # R_{i, y}^p

	crownArea = pi * crownRadius * crownRadius

	return (crownArea)
}

#### Define C0(p) and C1(p) table 2 Purves2007
C0_C1 = data.table(parameter = c("R0", "R40", "Rus", "B", "M", "Vus"),
	C0 = c(0.503, 0.5, 0.701, 0.196, 0.95, 2.551),
	C1 = c(3.126, 10, 3.955, 0.511, 0.95, 4.106))

getAllometries = function(species_id, allometries = purves2007_allometries)
{
	n = length(species_id)
	ls_allometries = list(a = numeric(length = n), b = numeric(length = n), T_param = numeric(length = n))
	for (i in 1:n)
	{
		ls_allometries$a[i] = allometries[species == species_id[i], a]
		ls_allometries$b[i] = allometries[species == species_id[i], b]
		ls_allometries$T_param[i] = allometries[species == species_id[i], T]
	}
	return (ls_allometries)
}

getAllometriesVec = function(species_id, allometries = purves2007_allometries)
{
	n = length(species_id)
	ls_species = unique(species_id)
	nb_sp = length(ls_species)

	ls_allometries = data.table(a = numeric(length = n), b = numeric(length = n), T_param = numeric(length = n))

	for (i in 1:nb_sp)
	{
		sp_ind = which(species_id == ls_species[i])
		ls_allometries[sp_ind, c("a", "b", "T_param") := .(allometries[species == ls_species[i], a],
			allometries[species == ls_species[i], b], allometries[species == ls_species[i], T])]
	}
	return (ls_allometries)
}

#### Check the crown radius of each individual at the estimated h*
crownRadius_star = function(height, distToTop, a, b, T_param, C0_C1)
{
	# Table S2
	R0_C0 = C0_C1[parameter == "R0", C0]
	R0_C1 = C0_C1[parameter == "R0", C1]

	R40_C0 = C0_C1[parameter == "R40", C0]
	R40_C1 = C0_C1[parameter == "R40", C1]

	M_C0 = C0_C1[parameter == "M", C0]
	M_C1 = C0_C1[parameter == "M", C1]

	B_C0 = C0_C1[parameter == "B", C0]
	B_C1 = C0_C1[parameter == "B", C1]

	# Appendix S3, Eq S3.3 (erroneously denoted S2.3 in the article)
	R0 = (1 - T_param)*R0_C0 + T_param*R0_C1
	R40 = (1 - T_param)*R40_C0 + T_param*R40_C1

	M = (1 - T_param)*M_C0 + T_param*M_C1
	B = (1 - T_param)*B_C0 + T_param*B_C1

	# All the following part is vectorised
	# Convert height to dbh. /!\ If height is in m, then dbh is in cm /!\
	dbh = heightToDbh(height, a, b)

	# Calculate potential max radius Eq S1.6
	Rp_max = R0 + (R40 - R0)*dbh/40

	n = length(height)
	crownRadius = Rp_max * ( pmin(distToTop, height*M) / (height*M) )^B # R_{i, y}^p

	return (crownRadius)
}
