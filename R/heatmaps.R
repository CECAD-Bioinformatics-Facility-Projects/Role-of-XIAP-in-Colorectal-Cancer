save_heatmap <- function(path, heatmap) {
	fs::dir_create(fs::path_dir(path))
	pdf(file = path, width = 10, height = 12)
	ComplexHeatmap::draw(heatmap)
	dev.off()
	return(path)
}
