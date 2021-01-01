


all: cpp root

cpp:
	@echo "Making C++ modules"
	make -C PulseQualityCheck_CppTemplate
	make -C LivetimeModule
	make -C InitializeRQFile 
	make -C PulseFinder_PODoverlap 
	make -C PulseTiming_HeightTiming
#	make -C Converters
	make -C PODSummer_PODoverlap
	make -C PulseQuantities_PseudoPromptFraction 
#	make -C PulseFinder_BlackBox
	make -C PulseCalibration_BaselineSanity
	make -C PulseCalibration_MakeCvtFile
	make -C PositionReconstruction_CorrCentroid
	make -C PositionReconstruction_Hitmap
	make -C MonteCarlo_AnswerKey
	make -C PulseQuantities_WaterPmtRQs 
	make -C PulseTiming_CompensationForTRC 
#	make -C PulseClassifier_PHET 

	-I /usr/common/usg/software/boost/1.54.0/include/
root:
	@echo "Making ROOT modules"
#	make -C PulseQuantities_FastMinimumSet

clean:
	@echo "Cleaning house!"
	make -C PulseQualityCheck_CppTemplate clean
	make -C LivetimeModule clean
	make -C InitializeRQFile clean
	make -C PulseFinder_PODoverlap clean
	make -C PulseTiming_HeightTiming clean
	make -C Converters clean
	make -C PODSummer_PODoverlap clean
	make -C PulseQuantities_PseudoPromptFraction clean
	make -C PulseFinder_BlackBox clean
	make -C PulseCalibration_BaselineSanity clean
	make -C PulseCalibration_MakeCvtFile clean
	make -C PulseQuantities_FastMinimumSet clean
	make -C PositionReconstruction_CorrCentroid clean
	make -C PulseQuantities_WaterPmtRQs clean
	make -C PulseTiming_CompensationForTRC clean
	make -C PulseClassifier_PHET clean
	make -C PositionReconstruction_Hitmap clean
	make -C MonteCarlo_AnswerKey clean

distclean:
	@echo "Some serious cleaning."
	rm -rf build/*.o
	rm -rf bin/*
