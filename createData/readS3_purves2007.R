
#### Aims of prog: to read the supplementary information from Purves 2007
#
## Remaks
# Rk1: Carya alba (tsn: 501306) = Carya tomentosa (tsn: 19247).
# Rk2: Chamaecyparis nootkatensis (tsn: 183451) = Callitropsis nootkatensis (tsn: 822596).
# Rk3: Quercus prinus (tsn: 19398) = Quercus montana (tsn: 19379)
# Rk4: Taxodium ascendens (tsn: 183433) = Taxodium distichum var. imbricarium (tsn: 541138)
# Rk5: Error in Table_S3.doc (Purves2007), Betula populifolla should be read Betula populifolia

#### Load packages
library(data.table)
library(textreadr)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()
options(max.print = 500)

#### Load data
appS3 = read_doc("../Purves2007_supp-meth/Table_S3.doc", skip = 0, remove.empty = TRUE, trim = TRUE, format = FALSE)
ls_species = readRDS("./ls_speciesParametrised.rds")
ls_species = ls_species[,1]

purves2007_allometries = readRDS("./purves2007_allometries.rds")

#### Tool function
## Sort species by alphabetical order
sortingSpecies = function(ls_species)
{
	tsn_dt = data.table(species = character(length(ls_species)),
		tsn = character(length(ls_species))) # tsn = Taxonomic Serial Number

	tsn_dt[, species := stri_sub(str = ls_species,
		from = stri_locate_first(str = ls_species, regex = "-")[,1] + 1)]

	tsn_dt[, tsn := stri_sub(str = ls_species, from = 1,
		to = stri_locate_first(str = ls_species, regex = "-")[,1] - 1)]
	tsn_dt = tsn_dt[order(species),]
	return (tsn_dt)
}

#### Vernacular/scientific names
## Data table tsn
tsn_dt = sortingSpecies(ls_species)
tsn_dt[, c("latin", "vernacular") := .(
	c("Abies amabilis", "Abies balsamea", "Abies concolor", "Abies grandis", "Abies lasiocarpa",
		"Abies magnifica", "Acer negundo", "Acer pensylvanicum", "Acer rubrum", "Acer saccharum",
		"Acer saccharinum", "Alnus rubra", "Betula alleghaniensis", "Betula lenta",
		"Betula papyrifera", "Betula populifolia", "Calocedrus decurrens", "Carya alba",
		"Carpinus caroliniana", "Carya cordiformis", "Carya glabra", "Carya ovata",
		"Carya texana", "Celtis laevigata", "Celtis occidentalis", "Cercocarpus ledifolius",
		"Chamaecyparis nootkatensis", "Cornus florida", "Diospyros virginiana", "Fagus grandifolia",
		"Fraxinus americana", "Fraxinus nigra", "Fraxinus pennsylvanica", "Ilex opaca",
		"Juglans nigra", "Juniperus ashei",
		"Juniperus occidentalis", "Juniperus osteosperma", "Juniperus scopulorum",
		"Juniperus virginiana", "Larix laricina", "Larix occidentalis", "Liquidambar styraciflua",
		"Liriodendron tulipifera", "Lithocarpus densiflorus", "Magnolia virginiana",
		"Nyssa aquatica", "Nyssa biflora", "Nyssa, sylvatica", "Ostrya virginiana",
		"Oxydendrum arboreum",
		"Picea engelmannii", "Picea glauca", "Picea mariana", "Picea rubens",
		"Picea sitchensis", "Pinus albicaulis", "Pinus banksiana", "Pinus contorta",
		"Pinus echinata", "Pinus edulis", "Pinus elliottii", "Pinus jeffreyi",
		"Pinus monticola", "Pinus palustris", "Pinus ponderosa", "Pinus resinosa",
		"Pinus strobus", "Pinus taeda", "Pinus virginiana", "Platanus occidentalis",
		"Populus balsamifera", "Populus grandidentata", "Populus tremuloides",
		"Prosopis glandulosa", "Prunus pensylvanica", "Prunus serotina",
		"Pseudotsuga menziesii", "Quercus alba", "Quercus chrysolepis", "Quercus coccinea",
		"Quercus ellipsoidalis", "Quercus falcata", "Quercus gambelii",
		"Quercus kelloggii", "Quercus laurifolia", "Quercus macrocarpa",
		"Quercus marilandica", "Quercus nigra", "Quercus phellos", "Quercus prinus",
		"Quercus rubra", "Quercus stellata", "Quercus velutina", "Quercus virginiana",
		"Robinia pseudoacacia", "Salix nigra", "Sassafras albidum", "Taxodium ascendens",
		"Taxodium distichum", "Thuja occidentalis", "Thuja plicata", "Tilia americana",
		"Tsuga canadensis", "Tsuga heterophylla", "Tsuga mertensiana",
		"Ulmus alata", "Ulmus americana", "Ulmus rubra"),
	c("Pacific silver fir", "Balsam fir", "White fir", "Grand fir", "Subalpine fir",
		"California red fir", "boxelder", "Stripped maple", "Red maple", "Sugar maple",
		"Silver maple", "Red alder", "Yellow birch", "Sweet birch",
		"White birch", "Grey birch", "Incense cedar", "Mockernut hickory",
		"American hornbeam", "Bitternut hickory", "Pignut hickory", "Shagbark hickory",
		"Black hickory", "Sugarberry", "Hackberry", "curlleaf cercocarpus",
		"Alaska cedar", "Flowering dogwood", "Common persimmon", "American beech",
		"White ash", "Black ash", "Green ash", "American holly",
		"Black wallnut", "Ashe juniper",
		"Western juniper", "Utah juniper", "Rocky mountain juniper",
		"Eastern redcedar", "Tamarack", "Western larch","Sweetgum",
		"Yellow-poplar", "Tanoak", "Sweetbay",
		"Water tupelo", "Swamp tupelo", "Blackgum", "Eastern hophombeam",
		"Sourwood",
		"Englemann spruce", "White spruce", "Black spruce", "Red spruce",
		"Sitka spruce", "Whitebark pine", "Jack pine", "Lodgepole pine",
		"Shortleaf pine", "Common pinyon", "Slash pine", "Jeffrey pine",
		"Western white pine", "Longleaf pine", "Ponderosa pine", "Red pine",
		"Eastern white pine", "Loblolly pine", "Virginiana pine", "Sycamore",
		"Balsam poplar", "Bigtooth aspen", "Quaking aspen",
		"Western honey mesquite", "Pin cherry", "Black cherry",
		"Douglas-fir", "White oak", "Canyon live oak", "Scarlet oak",
		"Northern pin oak", "Southern red oak", "Gambel oak",
		"California black oak", "Laurel oak", "Bur oak",
		"Blackjack oak", "Water oak", "Willow oak", "Chestnut oak",
		"Northern red oak", "Post oak", "Black oak", "Live oak",
		"Black locust", "Black willow", "Sassafras", "Pondcypress",
		"Baldcypress", "Northern white cedar", "Western redcedar", "Carolina basswood",
		"Eastern hemlock", "Western hemlock", "Mountain hemlock",
		"Winged elm", "American elm", "Slippery elm"))]

appS3 = stri_split_boundaries(appS3, type = "word", skip_word_none = TRUE, simplify = TRUE)

## List latin combinations
latin_matrix = stri_split_boundaries(tsn_dt[, latin], type = "word",
	skip_word_none = TRUE, simplify = TRUE)

ls_manualSp = integer(tsn_dt[, .N])

## Look for species in appendix S3
# tsn_dt[, c("a", "b", "t") := NULL]
tsn_dt[, c("a", "b", "t") := .(0, 0, 0)]

# Few species will through a warning => check them manually
for (i in 1:nrow(latin_matrix))
{
	ind_genus = which(appS3 == latin_matrix[i, 1])
	if (length(ind_genus) == 1)
		ind_params = c(2, 4, 5) + ind_genus

	if (length(ind_genus) > 1)
	{
		ind_species = which(appS3 == latin_matrix[i, 2])
		good_genus = which( (ind_genus + 1) %in% ind_species)

		if (length(good_genus) != 1)
		{
			ls_manualSp[i] = i
			print(paste0(i, ": ", latin_matrix[i, 1], " ", latin_matrix[i, 2]))
			next;
		}

		# if ( length(ind_genus) >= length(ind_species) )
			ind_params = c(2, 4, 5) + ind_genus[good_genus]

		# if ( length(ind_genus) < length(ind_species) )
			# ind_params = c(1, 3, 4) + ind_species[good_genus]
	}

	val = as.list(as.numeric(stri_split_boundaries(appS3[ind_params],
		type = "word", skip_word_none = TRUE, simplify = TRUE)))
	tsn_dt[i, c("a", "b", "t") := val]
}

ls_manualSp = ls_manualSp[ls_manualSp != 0]
purves2007_allometries[, c("a", "b", "T") := lapply(.SD, abs), .SDcols = c("a", "b", "T")]

setnames(x = purves2007_allometries, old = c("a", "b", "T"), new = c("a", "b", "t"))

# all.equal(tsn_dt[1:3, .(a, b, t)], purves2007_allometries[1:3, .(a, b, t)])

lsProb = unique(which(tsn_dt[, .(a, b, t)] != purves2007_allometries[, .(a, b, t)]) %% tsn_dt[, .N])

sub_tsn = tsn_dt[lsProb]
speciesProb = sub_tsn[a != 0, species]

purves2007_allometries[sp %in% speciesProb]
tsn_dt[species %in% speciesProb, .(species, a, b, t)]

#### Check manually the remaining species
purves2007_allometries[ls_manualSp]
