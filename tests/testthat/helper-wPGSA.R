data(GSE143371_count)
data(GSE143371_group)

GeneName1 <- c("FRG2C", "FRG2", "DUX4", "COL23A1", "MTRNR2L1")
GeneName2 <- c("FRG2C", "FRG2", "AATF")

dummy_GSE143371_count1 <- matrix(0, nrow=5, ncol=6)
dummy_GSE143371_count2 <- matrix(0, nrow=3, ncol=6)

rownames(dummy_GSE143371_count1) <- GeneName1
rownames(dummy_GSE143371_count2) <- GeneName2

SampleName <- c("Undiff_1", "Undiff_2", "Undiff_3", "Neu_1", "Neu_2", "Neu_3")
colnames(dummy_GSE143371_count1) <- SampleName
colnames(dummy_GSE143371_count2) <- SampleName

dummy_GSE143371_count1["FRG2C", ] <- c(1.406332, 1.558506, 1.346196, 1.236234, 1.229522, 1.238637)
dummy_GSE143371_count1["COL23A1", ] <- c(7.560481, 7.704230, 8.410127, 4.548904, 4.771929, 4.971121)
dummy_GSE143371_count1["MTRNR2L1", ] <- c(-1.470367, -1.505061, -1.503108, -1.470036, -1.506252, -1.498131)

dummy_GSE143371_count2["FRG2C", ] <- c(1.406332, 1.558506, 1.346196, 1.236234, 1.229522, 1.238637)
dummy_GSE143371_count2["AATF", ] <- c(9.583129, 9.607526, 9.795627, 9.579084, 9.696098, 9.580364)

dummy_GSE143371_count <- rbind(dummy_GSE143371_count1, dummy_GSE143371_count2)
dummy_network_data <- rbind(
		data.frame(
		TF_name="AATF_Others_BIU-87_SRX22582424",
		GeneName=GeneName1,
		weight=1),
		data.frame(
		TF_name="AATF_Blood_NALM-6_SRX2493169",
		GeneName=GeneName2,
		weight=1))
