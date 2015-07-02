/*
 * higher_order_SMC.cpp
 *
 *  Created on: 22-May-2015
 *      Author: nilxwam
 */

#include <unistd.h>
#include <iostream>
#include <string>
#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/log.h>
#include <barrett/standard_main_function.h>
//#include <boost/thread.hpp>
//#include <barrett/thread/null_mutex.h>

using namespace barrett;
using detail::waitForEnter;
#include <samlibs.h>

#include <Dynamics.hpp>
#include <SMC_higher_order.hpp>
#include <reference_signal.hpp>
#include <DMP.h>

template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm,
		systems::Wam<DOF>& wam) {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);


	Eigen::MatrixXd Coeff;
	Eigen::VectorXd Delta;
	Eigen::VectorXd Amp;
	Eigen::VectorXd Freq;
	Eigen::VectorXd StartPos;
	Eigen::MatrixXd A;
	Eigen::MatrixXd B;
	Eigen::VectorXd P;
	Eigen::VectorXd Q;

	Sam::initEigenMat<double>(A, Sam::readFile<double>("A.txt"));
	Sam::initEigenMat<double>(B, Sam::readFile<double>("B.txt"));
	Sam::initEigenVec<double>(P, Sam::readFile<double>("P.txt"));
	Sam::initEigenVec<double>(Q, Sam::readFile<double>("Q.txt"));

	Sam::initEigenMat<double>(Coeff, Sam::readFile<double>("coeff.txt"));
	Sam::initEigenVec<double>(Delta, Sam::readFile<double>("delta.txt"));
	Sam::initEigenVec<double>(Amp, Sam::readFile<double>("amp.txt"));
	Sam::initEigenVec<double>(Freq, Sam::readFile<double>("freq.txt"));
	Sam::initEigenVec<double>(StartPos, Sam::readFile<double>("start.txt"));

	typedef boost::tuple<double> tuple_type;
	typedef systems::TupleGrouper<double> tg_type;
	//typedef boost::tuple<double, jp_type, jv_type, ja_type, jp_type, jv_type, jt_type> tuple_type;
	//typedef systems::TupleGrouper<double, jp_type, jv_type, ja_type, jp_type,  jv_type, jt_type> tg_type;
	tg_type tg;
	char tmpFile[] = "btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}

	//tool_details<DOF> tool_values(wam);

	const double TRANSITION_DURATION = 0.5;
	double amplitude1, omega1;
	const Eigen::Matrix4d coeff = Coeff;
	const Eigen::Vector4d delta = Delta;
	jp_type startpos(0.0);


	startpos[0] = StartPos[0];
	startpos[1] = StartPos[1];
	startpos[2] = StartPos[2];
	startpos[3] = StartPos[3];

	float a[4] = { 25.0, 25.0, 25.0, 25.0 }; // PD values
	float b[4] = { 25.0 / 4.0, 25.0 / 4.0, 25.0 / 4.0, 25.0 / 4.0 }; // PD values
	float y0[4] = { -0.0231193,1.05807,0.00852946,2.20058}; // Initial state [x0,y0,z0]
	float goal[4] = { -0.0357198,1.03931,0.0283707,2.18311  };
	DMP<DOF> DMP(4, 500, a, b,  9.534, 0.05, y0, goal);
	DMP.CheckDMPGoal();
	cout<<"Goal = "<< DMP.IsGoal<<"\n";

	startpos[0] = y0[0];
	startpos[1] = y0[1];
	startpos[2] = y0[2];
	startpos[3] = y0[3];



	const Eigen::Vector4d JP_AMPLITUDE = Amp;
	const Eigen::Vector4d OMEGA = Freq;
	bool status = true;
	const Eigen::Matrix4d A1 = A;
	const Eigen::Matrix4d B1 = B;
	const float P1 = P[0];
	const float Q1 = Q[0];

	//J_ref<DOF> joint_ref(JP_AMPLITUDE, OMEGA, startpos);
	SMC_higher_order<DOF> slide(status, coeff, delta , A1, B1,P1, Q1 );
	Dynamics<DOF> nilu_dynamics;

	wam.gravityCompensate();
	printf("Press [Enter] to turn on torque control to go to zero position");
	waitForEnter();

	wam.moveTo(startpos);
	printf("Press [Enter] to turn on torque control to joint 2.");
	waitForEnter();
	printf("Error 1 \n");

	systems::Ramp time(pm.getExecutionManager(), 1.0);
	const size_t PERIOD_MULTIPLIER = 1;
	//systems::PeriodicDataLogger<tuple_type> logger(pm.getExecutionManager(),
	//		new log::RealTimeWriter<tuple_type>(tmpFile,
	//				PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
	//		PERIOD_MULTIPLIER);
	printf("Error 2 \n");

	//systems::connect(tg.output, logger.input);
	//systems::connect(time.output, joint_ref.timef);
	systems::connect(wam.jpOutput, slide.feedbackjpInput);
	systems::connect(wam.jvOutput, slide.feedbackjvInput);
	systems::connect(wam.jpOutput, nilu_dynamics.jpInputDynamics);
	systems::connect(wam.jvOutput, nilu_dynamics.jvInputDynamics);
	systems::connect(nilu_dynamics.MassMAtrixOutput, slide.M);
	systems::connect(nilu_dynamics.CVectorOutput, slide.C);
	//systems::connect(joint_ref.referencejpTrack, slide.referencejpInput);
	//systems::connect(joint_ref.referencejvTrack, slide.referencejvInput);
	//systems::connect(joint_ref.referencejaTrack, slide.referencejaInput);
	systems::connect(DMP.ref_jp, slide.referencejpInput);
	systems::connect(DMP.ref_jv, slide.referencejvInput);
	systems::connect(DMP.ref_ja, slide.referencejaInput);

	wam.trackReferenceSignal(slide.controlOutput);
	printf("Error 3 \n");

	//systems::connect(time.output, tg.template getInput<0>());
printf("Error 3.1 \n");
	//systems::connect(joint_ref.referencejpTrack, tg.template getInput<1>());
	//systems::connect(joint_ref.referencejvTrack, tg.template getInput<2>());
	//systems::connect(DMP.ref_jp, tg.template getInput<1>());
printf("Error 3.2 \n");
	//systems::connect(DMP.ref_jv, tg.template getInput<2>());
printf("Error 3.3 \n");
	//systems::connect(DMP.ref_ja, tg.template getInput<3>());
printf("Error 3.4 \n");
	//systems::connect(wam.jpOutput, tg.template getInput<4>());
printf("Error 3.5 \n");
	//systems::connect(wam.jvOutput, tg.template getInput<5>());
printf("Error 3.6 \n");
	//systems::connect(slide.controlOutput, tg.template getInput<6>());
	printf("Error 4 \n");
//systems::connect(wam.cpOutput, tg.template getInput<1>());
//systems::connect(wam.cvOutput, tg.template getInput<2>());
	time.smoothStart(TRANSITION_DURATION);
	printf("Press [Enter] to stop.");
	waitForEnter();
	printf("Error 4.1 \n");
cout<<"Goal = "<< DMP.IsGoal<<"\n";
printf("Error 4.1.1 \n");
	//disconnect(tg.template getInput<0>());printf("Error 4.1.1 \n");
	//disconnect(tg.template getInput<1>());printf("Error 4.1.2 \n");
	//disconnect(Tg.template getInput<2>());
	//disconnect(Tg.template getInput<3>());
	//disconnect(Tg.template getInput<4>());
	//disconnect(logger.input);printf("Error 4.1.3 \n");
	//std::cout<<logger.isLogging()<<"!\n";printf("Error 4.1.4 \n");
	//logger.closeLog();
	printf("Error 4.2 \n");
	//time.smoothStop(TRANSITION_DURATION);
	printf("Error 4.3 \n");
	disconnect(wam.input);printf("Error 4.4 \n");
	wam.idle();

	printf("Error 5 \n");

	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	//log::Reader<boost::tuple<tuple_type> > lr(tmpFile);
	//lr.exportCSV(argv[1]);
	printf("Error 6 \n");

	printf("Output written to %s.\n", argv[1]);
	//std::remove(tmpFile);

	return 0;
}
