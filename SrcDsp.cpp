// SrcDsp.cpp : Defines the entry point for the console application.
//

#ifdef _WIN32
#include "stdafx.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <ctime>
#include <iterator>
#include "coefficients.h"
#include <random>


#ifdef _WIN32
#include "..\..\SrcLibraries\SrcDsp\filters.h"
#include "..\..\SrcLibraries\SrcDsp\generators.h"
#include "..\..\SrcLibraries\SrcDsp\files.h"
#include "..\..\SrcLibraries\SrcDsp\upsampling_filters.h"
#include "..\..\SrcLibraries\SrcDsp\dnsampling_filters.h"
#include "..\..\SrcLibraries\SrcDsp\modulators.h"
#include "..\..\SrcLibraries\SrcDsp\correlators.h"
#include "..\..\SrcLibraries\SrcDsp\mixers.h"
#include "..\..\SrcLibraries\SrcDsp\demodulators.h"
#else
#include "../../src_libraries/dsp/filters.h"
#include "../../src_libraries/dsp/generators.h"
#include "../../src_libraries/dsp/files.h"
#include "../../src_libraries/dsp/upsampling_filters.h"
#include "../../src_libraries/dsp/dnsampling_filters.h"
#include "../../src_libraries/dsp/modulators.h"
#include "../../src_libraries/dsp/correlators.h"
#include "../../src_libraries/dsp/mixers.h"
#include "../../src_libraries/dsp/demodulators.h"
#endif

/*-----------------------------------------------------------------------------
Test SymbolMappingSdpsk<T>
@tparam T The output of the modulator is complex<T>

@param bitsFile File in which the bits used to modulate the waveform are stored
@param outFile File in which the output of the modulator is stored


------------------------------------------------------------------------------*/
template<class T>
bool testSymbolMapperSdpsk()
{
	using namespace std;
	using namespace dsptl;
	const size_t nbrBits = 100; // Number of bits to run through the modulator

	// Create the bit pattern
	vector<int8_t> bits(nbrBits);
	std::default_random_engine dre;
	std::uniform_int_distribution<int> di(0, 1);
	for (size_t index = 0; index < bits.size(); ++index)
		bits[index] = static_cast<int8_t>(di(dre));
	// Create the modulator object
	SymbolMapperSdpsk<T> mod{};
	// Create the ouput of the modulator
	vector<complex<T> > out(bits.size());
	// Run the modulator
	mod.step(bits, out);
	// Save the output file
	ofstream osBits("sdpsk_mapper_input.txt");
	ofstream osOut("sdsk_mapper_output.txt");
	saveAsciiSamples(bits, osBits);
	saveAsciiSamples(out, osOut);
	return false;
}


/*-----------------------------------------------------------------------------
Test SymbolMappingSdpsk<T>
@tparam T The output of the modulator is complex<T>

@param bitsFile File in which the bits used to modulate the waveform are stored
@param outFile File in which the output of the modulator is stored


------------------------------------------------------------------------------*/
template<class T>
bool testSymbolMapperQpsk()
{
	using namespace std;
	using namespace dsptl;
	const size_t nbrBits = 100; // Number of bits to run through the symbol mapper
	assert(nbrBits % 2 == 0);

	// Create the bit pattern
	vector<int8_t> bits(nbrBits);
	std::default_random_engine dre;
	std::uniform_int_distribution<int> di(0, 1);
	for (size_t index = 0; index < bits.size(); ++index)
		bits[index] = static_cast<int8_t>(di(dre));
	// Create the modulator object
	SymbolMapperQpsk<T> mapper{};
	// Create the ouput of the modulator
	vector<complex<T> > out(bits.size() / 2);
	// Run the modulator
	mapper.step(bits, out);
	// Save the output file
	ofstream osBits("qpsk_mapper_input.txt");
	ofstream osOut("qpsk_mapper_output.txt");
	saveAsciiSamples(bits, osBits);
	saveAsciiSamples(out, osOut);
	return false;
}


/*-----------------------------------------------------------------------------
Test the fixed pattern correlator

@arg The first test checks the dynamic range as well as the speed of the correlator
@arg The second test checks the correlator against canned waveform

------------------------------------------------------------------------------*/

template <class T>
bool testFixedPatternCorrelator()
{

	using namespace std;
	using namespace dsptl;
	const size_t nbrBits = 100; // Number of bits to run through the modulator

	//---------- Checks dynamic range as well as speed of the computation
	{
		// Create input vector
		vector<complex<T> > in(1000);
		for (int index = 0; index < static_cast<int>(in.size()); ++index)
			in[index] = complex<int32_t>(index, -index);  // Test for 14 bits quantizer

		// Create the correlator object
		FixedPatternCorrelator<>  corr{};

		//set the correlation pattern
		array<complex<int32_t>, 32> coeffs;
		for (size_t index = 0; index < coeffs.size(); ++index)
			coeffs[index] = complex<int32_t>(4000+ index+1, 0); // 14 bits correlation coefficients

		corr.setPattern(coeffs);

		// Run the correlator
		int corrIndex;


		cout << "\nStarting correlator speed test" << '\n';
		clock_t t = clock();
		int iterations = 100;
		for (int index = 0; index < iterations; index++)
			corr.step(in, corrIndex);
		t = clock() - t;
		cout << "correlator: Time per iteration: " << (((double)t) * 1000) / CLOCKS_PER_SEC / iterations << " milliseconds" << '\n';
	}
	//-----------------------------------------------------------------------
	//---------- Checks the correlator against a canned waveform
	//-----------------------------------------------------------------------
	{

		// Set the correlator coefficients
		// The file used to read the coefficients is a formatted file where each 
		// value is a double.
		// 0.61825828 -0.50000000
		// 0.87399734 - 0.61825828
		// 0.61825828 0.12092214
		// - 0.15550441 0.61825828

		const int corrSize = 32;
		string corrFile{ "DnlSyncOqpsk.dat" };
		ifstream iscorr{ corrFile };
		if (!iscorr)
		{
			cout << "\nCannot open correlator coefficient file " << corrFile << endl;
			exit(1);
		}
		vector<double> corrCoeffsDouble;
		copy(istream_iterator<double>(iscorr), istream_iterator<double >(), back_inserter(corrCoeffsDouble));
		assert(corrCoeffsDouble.size() == corrSize * 2);
		array<complex<int32_t>, corrSize> corrCoeffs;
		const int corrCoeffsScaleFactor = (INT16_MAX >> 3);
		for (int k = 0; k < corrSize; ++k)
		{
			int32_t a = static_cast<int32_t>(corrCoeffsDouble[2 * k] * corrCoeffsScaleFactor);
			int32_t b = static_cast<int32_t>(corrCoeffsDouble[2 * k + 1] * corrCoeffsScaleFactor);
			corrCoeffs[k] = { a, b };
		}

		// Create the correlator object

		FixedPatternCorrelator<>  corr{};
		corr.setPattern(corrCoeffs);


		// Prepare the input file
		// The file contains unformatted, interleaved signed 16 bits I,Q samples
		// The sampling rate is 38400 sps
		string inFile{ "nts_baseband_1.dat" };
		ifstream isinfile{ inFile, ios::in | ios::binary };
		if (!isinfile)
		{
			cout << "\nCannot open baseband samples file " << inFile << endl;
			exit(1);
		}

		// Read the file block by block and perform the correlation
		int corrIndex = 0;
		const int blockSize = 1024;
		vector<complex<int16_t>> inSamples(blockSize);
		char * ptr = reinterpret_cast<char *>(&inSamples[0]);
		while (isinfile.good())
		{
			isinfile.read(ptr, blockSize * sizeof(complex<int16_t>));
			// Note: the last block may have been filled only partially
			// so the last elements may be garbage. This is OK

			// For testing purposes, the inputs are scaled by 2 to see if we saturate somewhere
			int scale = 0;
			for (size_t index = 0; index < inSamples.size(); ++index)
			{
				inSamples[index] = complex<int16_t>(inSamples[index].real() >> scale, inSamples[index].imag() >> scale);

			}

			bool result = corr.step(inSamples, corrIndex);
			if (result)
			{
				cout << "\nCorrelation found at : " << corrIndex;
			}
			
		}


		// Save the output file
		//ofstream osBits(bitsFile);
		//ofstream osOut(outFile);
		//saveAsciiSamples(bits, osBits);
		//saveAsciiSamples(out, osOut);
	}
	return false;
}

template <class InType, class OutType, class InternalType, class CoefType = int32_t>
bool testFilters()
{
	using namespace std;
	using namespace dsptl;

	// Impulse response test
	{
		// Impulse response test
		cout << "Creation of filter object" << '\n';
		vector<CoefType> filterCoeff{ 134 , 309, 768};
		FilterFir<InType, OutType, InternalType, CoefType> filter(filterCoeff);

		cout << "Filter inpulse response" << '\n';
		vector<InType> input{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		vector<OutType> output(input.size());
		filter.step(input, output);

		for (size_t index = 0; index < output.size(); ++index)
			cout << output[index] << " ";

		cout << '\n';

		cout << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);


		for (size_t index = 0; index < output.size(); ++index)
			cout << output[index] << " ";
		cout << '\n';
	}

	// Sinewave Test
	{
		// Create input and output vectors
		int outLength = 500;
		vector<InType> in(outLength);
		vector<OutType> out(outLength);



		// Create the vector for the coefficients
		vector<int32_t> coeffs;
		for (size_t index = 0; index < coeffsRrc.size(); ++index)
		{
			coeffs.push_back(static_cast<int32_t>(coeffsRrc[index] * (INT16_MAX >> 0)));
		}
		// Create filter object
		FilterFir<InType, OutType, InternalType, CoefType> filter(coeffs);

		// Create debugging vectors and files 

		string sineFilename{ "filter_input.txt" };
		string filterOutFilename{ "filter_output.txt" };
		ofstream osIn(sineFilename);
		ofstream osOut(filterOutFilename);
		// vectors to store the results of all iterations
		vector<InType> inCombined;
		vector<OutType> outCombined;
		filter.reset();
		GenSine<InType> gen(0.1, 32000);
		int nbrLoops = 5;
		for (int loopIndex = 0; loopIndex < nbrLoops; ++loopIndex)
		{
			gen.step(in);
			filter.step(in, out);
			inCombined.insert(inCombined.end(), in.cbegin(), in.cend());
			outCombined.insert(outCombined.end(), out.cbegin(), out.cend());

		}

		saveAsciiSamples(outCombined, osOut);
		saveAsciiSamples(inCombined, osIn);
		cout << "Filter input saved in :" << sineFilename << '\n';
		cout << "Filter output saved in : " << filterOutFilename << '\n';

	}



	return false;
}


template < class InType, class OutType, class PhaseType, unsigned N = 4096 >
bool testMixers()
{
	using namespace std;
	using namespace dsptl;

	// Sinewave Test
	{
		// Create input and output vectors
		int outLength = 500;
		vector<InType> in(outLength);
		vector<OutType> out(outLength);

		// Create mixer object
		Mixer< InType,OutType,PhaseType, N > mixer;

		// Create debugging vectors and files 

		string sineFilename{ "mixer_input.txt" };
		string filterOutFilename{ "mixer_output.txt" };
		ofstream osIn(sineFilename);
		ofstream osOut(filterOutFilename);
		// vectors to store the results of all iterations
		vector<InType> inCombined;
		vector<OutType> outCombined;
		mixer.reset();
		mixer.setFrequency(0.9f);
		GenSine<InType> gen(0.05, 32000);
		int nbrLoops = 5;
		for (int loopIndex = 0; loopIndex < nbrLoops; ++loopIndex)
		{
			gen.step(in);
			mixer.step(in, out);
			inCombined.insert(inCombined.end(), in.cbegin(), in.cend());
			outCombined.insert(outCombined.end(), out.cbegin(), out.cend());

		}

		saveAsciiSamples(outCombined, osOut);
		saveAsciiSamples(inCombined, osIn);
		cout << "Mixer input saved in :" << sineFilename << '\n';
		cout << "Mixer output saved in : " << filterOutFilename << '\n';

	}



	return false;
}


/*-----------------------------------------------------------------------------
Test a downsampling filter.

The template parameters mirrors the template parameters available for the 
filter.

@tparam InType Type of the input signal. Can be float, double, complex, int...
@tparam OutType Type of the output signal
@tparam InternalType Type used internally for the computation
@tparam CoefType Type of the coefficients. Cannot be complex
@tparam M	Downsampling ratio

It is the responsibility of the caller to make sure that the different types
work smoothly. Overflow and underflow conditions must not occur

@return false if no error. true if an error occurred

------------------------------------------------------------------------------*/
template<class InType, class OutType, class InternalType, class CoefType, unsigned M>
bool testDnsamplingFilter()
{
	using namespace dsptl;
	using namespace std;

	vector<InType> in;
	vector<OutType> out;


	// The following coefficients correspond to:
	// Lowpass - Fpass = 0.15 - Fstop = 0.4 - Apass = 0.5 - Astop = 50
	const vector<double> coeffs_double = {
		-0.01059321924411, -0.02235300003302, -0.02449850003147, 0.001614827776104,
		0.06449213811506, 0.1507777255143, 0.2269893874909, 0.2574812935042,
		0.2269893874909, 0.1507777255143, 0.06449213811506, 0.001614827776104,
		-0.02449850003147, -0.02235300003302, -0.01059321924411
	};



	vector<CoefType> coeffs;
	for (auto it = coeffs_double.cbegin(); it != coeffs_double.cend(); ++it)
		coeffs.push_back(static_cast<CoefType>((*it) * (INT16_MAX >> 0)));


	FilterDnsamplingFir<InType, OutType, InternalType, CoefType, M> filter(coeffs);
	size_t outLength ;


	// Impulse response test
	// The output should be a decimated version of the coefficients
	outLength = 500;
	in.assign(M * outLength, {});
	in[0] = { 1, 1 };
	out.assign(outLength, {});
	filter.reset();
	filter.step(in, out);

	// History storage test
	// The internal history buffer should contain the
	// last N-1 values of the input.
	outLength = 500;
	in.assign(M * outLength, {});
	for (size_t k = 0; k < in.size(); ++k)
		in[k] = InType(k, k);
	out.assign(outLength, {});
	filter.reset();
	filter.step(in, out);

	// Sinewave test
	// This test calls the step function of the filter
	// several times to make sure that the filter history is
	// handled properly
	outLength = 500;
	in.assign(M * outLength, {});
	out.assign(outLength, {});
	filter.reset();


	string sineFilename{ "filter_dn_input.txt" };
	string filterOutFilename{ "filter_dn_output.txt" };
	ofstream osIn(sineFilename);
	ofstream osOut(filterOutFilename);
	// vectors to store the results of all iterations
	vector<InType> inCombined;
	vector<OutType> outCombined;

	GenSine<InType> gen(0.1, 150);
	int nbrLoops = 5;
	for (int loopIndex = 0; loopIndex < nbrLoops; ++loopIndex)
	{
		gen.step(in);
		filter.step(in, out);
		inCombined.insert(inCombined.end(), in.cbegin(), in.cend());
		outCombined.insert(outCombined.end(), out.cbegin(), out.cend());

	}

	saveAsciiSamples(outCombined, osOut);
	saveAsciiSamples(inCombined, osIn);
	cout << "Filter input saved in :" << sineFilename << '\n';
	cout << "Filter output saved in : " << filterOutFilename << '\n';



	return false;
}


/*-----------------------------------------------------------------------------
Test an upsampling filter

The template parameters mirrors the template parameters available for the 
filter.\n
The routine supports decimation by 8 and by 2. Support to other upsampling ratio
can be added by adding additional filters.

@tparam InType Type of the input signal. Can be float, double, complex, int...
@tparam OutType Type of the output signal
@tparam InternalType Type used internally for the computation
@tparam CoefType Type of the coefficients. Cannot be complex
@tparam L	Upsampling ratio

It is the responsibility of the caller to make sure that the different types
work smoothly. Overflow and underflow conditions must not occur

@return false if no error. true if an error occurred

------------------------------------------------------------------------------*/
template<class InType, class OutType, class InternalType, class CoefType, unsigned L = 2>
bool testUpsamplingFilters()
{
	using namespace std;
	using namespace dsptl;

	assert(L == 2 || L == 8);

	cout << "\n ********* Upsampling Filters Test  ******** " << "\n\n";

	// Filter with a pass band from 0 to 0.4 and a stop band from 0.5 to 1
	// Linear phase
	// Passband ripple is 0.5 dB
	// Stopband attenuation is 80 dB
	// The low frequency gain is nominally 0
	// Filter was designed with the equiripple method
	// This filter can be used for upsampling by 2
	// One zero is added to make it a multiple of 2
	
	const vector<double> filterCoeffUp2{
		   1.02715039455e-05,-0.0009507952728189,-0.003197553296403,-0.004799208621488,
		  -0.002780826822133, 0.002022711644296, 0.003872174862249,-0.0004787585196786,
		  -0.005137122479716,-0.001982060672696, 0.005736105386769,  0.00549530744134,
		  -0.004789088290104,-0.009562021005043, 0.001580427155189,  0.01322572211686,
		   0.004292340148608, -0.01516534108286, -0.01278102298903,  0.01380490128625,
			0.02332688888775,-0.007341801410271, -0.03489034785084,-0.006549609457005,
			0.04611214859014,  0.03238686628952, -0.05550459317353, -0.08581806818547,
			0.06174961485047,    0.311280608227,   0.4360604051562,    0.311280608227,
			0.06174961485047, -0.08581806818547, -0.05550459317353,  0.03238686628952,
			0.04611214859014,-0.006549609457005, -0.03489034785084,-0.007341801410271,
			0.02332688888775,  0.01380490128625, -0.01278102298903, -0.01516534108286,
		   0.004292340148608,  0.01322572211686, 0.001580427155189,-0.009562021005043,
		  -0.004789088290104,  0.00549530744134, 0.005736105386769,-0.001982060672696,
		  -0.005137122479716,-0.0004787585196786, 0.003872174862249, 0.002022711644296,
		  -0.002780826822133,-0.004799208621488,-0.003197553296403,-0.0009507952728189,
		   1.02715039455e-05,0};

	// Filter with a pass band from 0 to 0.100 and a stop band from 0.125 to 1
	// Linear phase
	// Passband ripple is 0.5 dB
	// Stopband attenuation is 80 dB
	// The low frequency gain is nominally 0
	// Filter was designed with the equiripple method
	// This filter can be used for upsampling by 8
	// 6 zeros are added to make it a multiple of 8


	const vector<double> filterCoeffUp8{
		  5.833238323783e-05,-1.557760923617e-05,-4.719868618434e-05,-0.0001020146992172,
		  -0.0001817133043676,-0.0002865581821126,-0.0004144512039692,-0.0005603576734241,
		  -0.0007161033590357,-0.0008705731999149,-0.001010327623762,-0.001120650661823,
		  -0.001186982240667,-0.001196586407912,-0.001140252017006,-0.001013822640645,
		  -0.0008194041527966,-0.0005661025859619,-0.0002700266216707,4.667139003715e-05,
		   0.000357535562294,0.0006336393775786,0.0008472142381677,0.0009741576180476,
		  0.0009973930490773,0.0009093102706575,0.0007135188617084, 0.000425480020635,
		  7.184197604252e-05,-0.0003115902776905,-0.0006831301254212,-0.0009994524953701,
		  -0.001220398163793, -0.00131383744119,-0.001260031144097,-0.001054887783767,
		  -0.0007116204689155,-0.0002605415630515,  0.00025314023239,0.0007735395181854,
		   0.001239924697254, 0.001593379701041, 0.001784098020321, 0.001777820218265,
		   0.001561146301753, 0.001144649819349,0.0005633939642417,-0.0001254511161995,
		  -0.0008478025039868,-0.001520268110295,-0.002059425158352,-0.002391786398817,
		   -0.00246329210467,-0.002247146520675,-0.001748969262069,-0.001008439785949,
		  -9.69455133972e-05,0.0008887313819455,  0.00183618137967, 0.002629734350894,
		   0.003164099532977, 0.003357618776919, 0.003163756658381, 0.002579040707362,
		   0.001646484215544,0.0004536448083273,-0.0008747658910921,-0.002189230416327,
		  -0.003331585495313,-0.004153325160225,-0.004533998026211, -0.00439752935277,
		  -0.003724378402722,-0.002557777209232,-0.001002929933222,0.0007811982239899,
		   0.002597066458482, 0.004229377436137, 0.005469419495592,  0.00614034751289,
		   0.006120472853179, 0.005361615090755, 0.003900157956269,  0.00185886619426,
		  -0.0005612863484069,-0.003099066019874,-0.005459280602421,-0.007344964969815,
		  -0.008492144805135,-0.008703360641772,-0.007875991457107,-0.006021783113354,
		  -0.003274736719285,0.0001144012537817, 0.003797383848055, 0.007360268155712,
			0.01036518538984,  0.01239818639022,  0.01311833901868,  0.01230288715816,
		   0.009883343186466, 0.005968000267456,0.0008476552668364,-0.005017234439273,
		   -0.01102778663147, -0.01649428357353, -0.02069629015272, -0.02294997751575,
		   -0.02267636876049, -0.01946385882907, -0.01311860279437,-0.003697348553269,
		   0.008481178646615,  0.02284831106083,  0.03862062098673,  0.05485571572066,
			0.07052410627539,  0.08459094656979,  0.09610034117727,   0.1042546073867,
			  0.108481221273,    0.108481221273,   0.1042546073867,  0.09610034117727,
			0.08459094656979,  0.07052410627539,  0.05485571572066,  0.03862062098673,
			0.02284831106083, 0.008481178646615,-0.003697348553269, -0.01311860279437,
		   -0.01946385882907, -0.02267636876049, -0.02294997751575, -0.02069629015272,
		   -0.01649428357353, -0.01102778663147,-0.005017234439273,0.0008476552668364,
		   0.005968000267456, 0.009883343186466,  0.01230288715816,  0.01311833901868,
			0.01239818639022,  0.01036518538984, 0.007360268155712, 0.003797383848055,
		  0.0001144012537817,-0.003274736719285,-0.006021783113354,-0.007875991457107,
		  -0.008703360641772,-0.008492144805135,-0.007344964969815,-0.005459280602421,
		  -0.003099066019874,-0.0005612863484069,  0.00185886619426, 0.003900157956269,
		   0.005361615090755, 0.006120472853179,  0.00614034751289, 0.005469419495592,
		   0.004229377436137, 0.002597066458482,0.0007811982239899,-0.001002929933222,
		  -0.002557777209232,-0.003724378402722, -0.00439752935277,-0.004533998026211,
		  -0.004153325160225,-0.003331585495313,-0.002189230416327,-0.0008747658910921,
		  0.0004536448083273, 0.001646484215544, 0.002579040707362, 0.003163756658381,
		   0.003357618776919, 0.003164099532977, 0.002629734350894,  0.00183618137967,
		  0.0008887313819455,-9.69455133972e-05,-0.001008439785949,-0.001748969262069,
		  -0.002247146520675, -0.00246329210467,-0.002391786398817,-0.002059425158352,
		  -0.001520268110295,-0.0008478025039868,-0.0001254511161995,0.0005633939642417,
		   0.001144649819349, 0.001561146301753, 0.001777820218265, 0.001784098020321,
		   0.001593379701041, 0.001239924697254,0.0007735395181854,  0.00025314023239,
		  -0.0002605415630515,-0.0007116204689155,-0.001054887783767,-0.001260031144097,
		   -0.00131383744119,-0.001220398163793,-0.0009994524953701,-0.0006831301254212,
		  -0.0003115902776905,7.184197604252e-05, 0.000425480020635,0.0007135188617084,
		  0.0009093102706575,0.0009973930490773,0.0009741576180476,0.0008472142381677,
		  0.0006336393775786, 0.000357535562294,4.667139003715e-05,-0.0002700266216707,
		  -0.0005661025859619,-0.0008194041527966,-0.001013822640645,-0.001140252017006,
		  -0.001196586407912,-0.001186982240667,-0.001120650661823,-0.001010327623762,
		  -0.0008705731999149,-0.0007161033590357,-0.0005603576734241,-0.0004144512039692,
		  -0.0002865581821126,-0.0001817133043676,-0.0001020146992172,-4.719868618434e-05,
		  -1.557760923617e-05,5.833238323783e-05,0,0,0,0,0,0
		};

	// ---- Filter coefficient selection
	vector<double> filterCoeffDouble;
	if (L == 2) filterCoeffDouble = filterCoeffUp2;
	else if (L == 8) filterCoeffDouble = filterCoeffUp8;
	assert(filterCoeffDouble.size() % L == 0);
	vector<CoefType> filterCoeff;
	for (auto it = filterCoeffDouble.cbegin(); it != filterCoeffDouble.cend(); ++it)
	{
		// filter coefficients are scaled by the size of the input type
		filterCoeff.push_back(static_cast<CoefType>(*it * numeric_limits<typename InType::value_type>::max()));
	}


	// ---- Display of filter parameters
	cout << "Upsampling Ratio: " << L << '\n';
	cout << "Number of coefficients: " << filterCoeff.size() << '\n';



	// ----- Sinewave response
	int samplesPerLoop = 100;
	cout << "\nSinewave response" << '\n';
	string inFilename{ "filter_up_input.txt" };
	string outFilename{ "filter_up_output.txt" };
	ofstream osIn(inFilename);
	ofstream osOut(outFilename);
	// vectors to store the results of all iterations
	vector<InType> inCombined;
	vector<OutType> outCombined;
	// Iteration vectors 
	vector<InType > in(samplesPerLoop);
	vector<OutType > out;
	out.resize(in.size() * L);
	// Objects creation
	GenSine<InType> gen(0.1, numeric_limits<typename InType::value_type>::max());
	FilterUpsamplingFir<InType, OutType, InternalType, CoefType, L> filter(filterCoeff);
	filter.reset();

	int nbrLoops = 2;
	// Processing loop
	for (int index = 0; index < nbrLoops; ++index)
	{
		gen.step(in);
		filter.step(in, out);
		inCombined.insert(inCombined.end(), in.cbegin(), in.cend());
		outCombined.insert(outCombined.end(), out.cbegin(), out.cend());
	}

	// ---- Save the results in a file
	saveAsciiSamples(outCombined, osOut);
	saveAsciiSamples(inCombined, osIn);
	cout << "Filter input saved in :" << inFilename << '\n';
	cout << "Filter output saved in : " << outFilename << '\n';


	return false;
}


/*-----------------------------------------------------------------------------
This routine test the OQPSKdemodulator by comparing the behavior of the routine
versus its matlab representation. \n

Three files are used :
@arg DemodulatorRefData: bit samples of the sync word with modulation removed.
@arg DemodulatorInputData.dat: Input of the demodulator exactly as Matlab sees it
@arg DemodulatorSoftbits.dat: Softbits output of the Matlab demdulator for a complete frame

As the rounding in Matlab may be slightly different than in C++, the softbbits may differ by
1 or 2. In addtion, Matlab softbits may exceed 127 and -128 while the C++ softbits saturate at 127 and -128\n

All the matlab data was created with an initial frequency offset of 0\n

The resulting error count should be zero

@return false if no error. true if an error occurred

------------------------------------------------------------------------------*/
bool testDemodulatorOqpsk()
{
	using namespace std;
	cout << "DemodulatorOqpsk Test" << '\n';
	const int syncSize = 32;

	// Create demodulator object
	dsptl::DemodulatorOqpsk<int16_t> demod;
	// Bit sync pattern
	vector<int8_t> bitSyncPattern = { 1,0,1,1,0,0,1,1,1,0,0,1,1,0,0,0,1,0,0,1,1,0,1,1,0,0,1,0,0,0,1,1 };
	assert(bitSyncPattern.size() == 32);
	demod.setSyncPattern(bitSyncPattern);
	// Initial frequency and phase of the loop
	demod.setInitialFrequency(0);
	demod.setInitialPhase(0);


	// -----  Read the reference vector from a file

	// Read the samples corresponding to the reference from  a file
	// the format is formatted I Q <RC> where I and Q are integers
	string refFile{ "DemodulatorRefData.txt" };
	ifstream isref{ refFile };
	if (!isref)
	{
		cout << "\nCannot open demodulator reference file " << refFile << endl;
		exit(1);
	}

	vector<int16_t> tmpVector;
	// Copy first the whole file content as a sequence of integers
	copy(istream_iterator<int16_t>(isref), istream_iterator<int16_t >(), back_inserter(tmpVector));
	assert(tmpVector.size() == syncSize * 2);
	// Reorganize the data as a vector of complex
	vector<complex<int16_t> > refVector;
	const int refVectorScaleFactor = 0;
	for (int k = 0; k < syncSize ; ++k)
	{
		refVector.push_back(complex<int16_t>(tmpVector[2 * k] >> refVectorScaleFactor, tmpVector[2 * k + 1] >> refVectorScaleFactor));
	}
	assert(refVector.size() == syncSize);
	demod.setReference(refVector);
	demod.reset();

	// -----  Prepare the input file

	// Prepare the input file
	// The file contains unformatted, interleaved signed 16 bits I,Q samples
	// there is one sample per bit.
	string inFile{ "DemodulatorInputData.dat" };
	ifstream isinfile{ inFile, ios::in | ios::binary };
	if (!isinfile)
	{
		cout << "\nCannot open demodulator input file " << inFile << endl;
		exit(1);
	}


	// -----  Perform iteration through the whole file content


	// Read the file block by block and perform the demodulation
	vector<int8_t> allBits;
	const int blockSize = 256;
	vector<complex<int16_t>> inSamples(blockSize);
	char * ptr = reinterpret_cast<char *>(&inSamples[0]);
	int32_t error;
	while (isinfile.good())
	{
		isinfile.read(ptr, blockSize * sizeof(complex<int16_t>));
		// Note: the last block may have been filled only partially
		// so the last elements may be garbage. This is OK
		vector<int8_t> out = demod.step(inSamples, error);
		// Append the temporary results to the previous results
		allBits.insert(allBits.end(), out.begin(), out.end());

	}


	// -----  Compare the result with the expected result

	// Read the softbits resulting from the result of the processing of
	// the first frame of the matlab demodulator
	string softbitsFile{ "DemodulatorOutputSoftbits1.txt" };
	ifstream issoftbits{ softbitsFile };
	if (!issoftbits)
	{
		cout << "\nCannot open softbits file " << softbitsFile << endl;
		exit(1);
	}

	vector<int16_t> matlabSoftbits;
	// Copy first the whole file content as a sequence of integers
	copy(istream_iterator<int16_t>(issoftbits), istream_iterator<int16_t >(), back_inserter(matlabSoftbits));
	assert(matlabSoftbits.size() == 9600 - 32);

	// For each softbits in the matlab file, we compare with the softbit
	// recovered by the C code.
	assert(allBits.size() >= matlabSoftbits.size());
	int errorCount = 0;
	for (size_t index = 0; index < matlabSoftbits.size(); ++index)
	{
		if (abs(matlabSoftbits[index] - allBits[index]) > 3)
		{
			// Potential for error. However we check if the matlab
			// was exceeding the 8 bits value
			if (matlabSoftbits[index] > 127 && allBits[index] == 127)
				continue;
			else if (matlabSoftbits[index] < -128 && allBits[index] == -128)
				continue;
			++errorCount;
			cout << "Matlab: " << matlabSoftbits[index] << "  C++ : " << static_cast<int16_t>(allBits[index]) << '\n';
		}
	}
	cout << "Error Count: " << errorCount << '\n';




	return errorCount == 0;
}


bool testFilters();
bool testGenerators();
bool testFiles();
bool testUpsamplingFilters();
int common_main();

#ifdef _WIN32
int _tmain(int argc, _TCHAR* argv[])
{
	return common_main();
}
#else
int main(int argc, char ** argv)
{
	return common_main();
}

#endif

int common_main()
{
	using namespace dsptl;
	using namespace std;

//	testDnsamplingFilter<complex<int16_t>, complex<int16_t>, complex<int32_t>, int32_t, 2>();
	//testSymbolMapperQpsk<int16_t>();

//	testFixedPatternCorrelator<int16_t>();

//	testFiles();
//	testGenerators();
	//testDemodulatorOqpsk();
	//testFilters<complex<int16_t>, complex<int16_t>, complex<int32_t>, int32_t>();
	//testMixers<std::complex<int16_t>, std::complex<int16_t>, int16_t, 4096 >();
	testUpsamplingFilters<complex<int16_t>, complex<int16_t>, complex<int32_t>, int32_t, 8>();

//	testSymbolMapperSdpsk<float>();
//	testSymbolMapperSdpsk<double>();
	testSymbolMapperSdpsk<int16_t>();
//	testSymbolMapperSdpsk<int32_t>();

	return 0;
}


bool testFiles()
{
	using namespace dsptl;
	using namespace std;

	{
		std::vector<int16_t> in{ 12, 34, 56, 78, 1, 43 };
		ofstream os("file_int16.txt");
		saveAsciiSamples(in, os);
		os.close();
	}

	{
		std::vector<float> in{ 12.34f, 34.5f, 56.67f, 78.f, 1.45f, 43.f };
		ofstream os("file_float.txt");
		saveAsciiSamples(in, os);
		os.close();
	}

	{
		std::vector<int16_t> in{ 12, 34, 56, 78, 1, 43 };
		ofstream os("file_int16.bin");
		saveBinarySamples(in, os);
		os.close();
	}


	{
		std::vector<complex<int16_t> > in{ { 12, 2 }, { 34, 56 }, { 78, 1 }, { 43, -14 } };
		ofstream os("file_complexint16.txt");
		saveAsciiSamples(in, os);
		os.close();
	}

	return false;
	

}


bool testGenerators()
{
	using namespace dsptl;
	using namespace std;

	{
#ifdef CPLUPLUS11
		using valueType = double;
#else
		typedef double valueType;
#endif
		GenSine<valueType> gen(0.4, 34.5);
		vector<valueType> out(20);

		gen.step(out);

#ifdef CPLUSPLUS11
		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
#else
		for(size_t index = 0; index < out.size(); index++)
			cout << out[index] << ' ';
		cout << "\n\n";
#endif
	}

	{
		typedef int16_t valueType ;
		GenSine<valueType> gen(1, 1000);
		vector<valueType> out(20);

		gen.step(out);

#ifdef CPLUSPLUS11
		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
#else
		for(size_t index = 0; index < out.size(); index++)
			cout << out[index] << ' ';
		cout << "\n\n";
#endif
	}


	{
#ifdef CPLUSPLUS11
		using valueType = complex<int16_t>;
#else
		typedef complex<int16_t> valueType;
#endif
		GenSine<valueType> gen(1, 1000);
		vector<valueType> out(20);

		gen.step(out);

#ifdef CPLUSPLUS11
		for (const auto& elt : out)
			cout << elt << ' ';
		cout << "\n\n";
#else
		for(size_t index = 0; index < out.size(); index++)
			cout << out[index] << ' ';
		cout << "\n\n";
#endif
	}


	return false;

}

#if 0

bool testUpsamplingFilters()
{
	using namespace std;
	using namespace dsptl;

	cout << "\n ********* Upsampling Filters Test  ******** " << "\n\n";
	const unsigned ratio = 2; // Upsampling ratio

	{
		cout << "--- Real Floating point test" << '\n';

		cout << "Creation of filter object" << '\n';
		// Filter with a transition baqndwidth between 0.25 and 0.5 normalized frequency at highest rate
		// 1 is Fs/2 or pi rad/samples
		vector<double> filterCoeff{
			0.009821820708929, -0.04103122121551, -0.05751222300691, 0.01999189052168,
			0.1968631695768, 0.3525450588422, 0.3525450588422, 0.1968631695768,
			0.01999189052168, -0.05751222300691, -0.04103122121551, 0.009821820708929
		};


		FilterUpsamplingFir<float, double, double, double, ratio> filter(filterCoeff);

		// ---- Display of filter parameters
		cout << "Upsampling Ratio: " << ratio << '\n';
		for (size_t index = 0; index< filterCoeff.size(); ++index)
			cout << filterCoeff[index] << " ";
		cout << '\n';

		// ----- Impulse response
		cout << "Filter inpulse response" << '\n';
		vector<float> input{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		vector<double> output(input.size() * ratio);
		filter.reset();
		filter.step(input, output);
		cout << "Input:\n";
		for (size_t index = 0; index< input.size(); ++index)
			cout << input[index] << " ";
		cout << '\n';
		cout << "Output:\n";
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
		cout << '\n';

		// ----- Delayed impulse response
		cout << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size() * ratio);
		filter.reset();
		filter.step(input, output);
		cout << "Input:\n";
		for (size_t index = 0; index< input.size(); ++index)
			cout << input[index] << " ";
		cout << '\n';
		cout << "Output:\n";
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
		cout << '\n';

		// ----- Sinewave response
		cout << "Sinewave response" << '\n';
		GenSine<float> gen(2.06, 34.5);
		vector<float> sineSignal(1000);
		gen.step(sineSignal);
		output.resize(sineSignal.size() * ratio);
		filter.reset();
		filter.step(sineSignal, output);
		string sineFilename{ "sine_input.txt" };
		string filterOutFilename{ "filter_output.txt"};
		ofstream osIn(sineFilename);
		ofstream osOut(filterOutFilename);
		saveAsciiSamples(output, osOut);
		saveAsciiSamples(sineSignal, osIn);
		cout << "Filter input saved in :" << sineFilename << '\n';
		cout << "Filter output saved in : " << filterOutFilename << '\n';

	}

	{
		cout << "--- Complex Double Floating point test" << '\n';

		cout << "Creation of filter object" << '\n';
		// Filter with a transition baqndwidth between 0.25 and 0.5 normalized frequency at highest rate
		// 1 is Fs/2 or pi rad/samples
		vector<double> filterCoeff{
			0.009821820708929, -0.04103122121551, -0.05751222300691, 0.01999189052168,
			0.1968631695768, 0.3525450588422, 0.3525450588422, 0.1968631695768,
			0.01999189052168, -0.05751222300691, -0.04103122121551, 0.009821820708929
		};


		FilterUpsamplingFir<complex<double>, complex<double>, complex<double>, double, ratio> filter(filterCoeff);

		// ---- Display of filter parameters
		cout << "Upsampling Ratio: " << ratio << '\n';
		for (size_t index = 0; index< filterCoeff.size(); ++index)
			cout << filterCoeff[index] << " ";
		cout << '\n';

		// ----- Impulse response
		cout << "Filter inpulse response" << '\n';
		vector<complex<double> > input{ { 1, 1 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };
		vector<complex<double> > output(input.size() * ratio);
		filter.reset();
		filter.step(input, output);
		cout << "Input:\n";
		for (size_t index = 0; index< input.size(); ++index)
			cout << input[index] << " ";
		cout << '\n';
		cout << "Output:\n";
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
		cout << '\n';


		// ----- Sinewave response
		cout << "Sinewave response" << '\n';
		GenSine<complex<double> > gen(2.06, 34.5);
		vector<complex<double> > sineSignal(1000);
		gen.step(sineSignal);
		output.resize(sineSignal.size() * ratio);
		filter.reset();
		filter.step(sineSignal, output);
		string sineFilename{ "sine_input.txt" };
		string filterOutFilename{ "filter_output.txt" };
		ofstream osIn(sineFilename);
		ofstream osOut(filterOutFilename);
		saveAsciiSamples(output, osOut);
		saveAsciiSamples(sineSignal, osIn);
		cout << "Filter input saved in :" << sineFilename << '\n';
		cout << "Filter output saved in : " << filterOutFilename << '\n';

	}

	{
		cout << "--- Complex Integer 32 test" << '\n';

		cout << "Creation of filter object" << '\n';
		// Filter with a transition baqndwidth between 0.25 and 0.5 normalized frequency at highest rate
		// 1 is Fs/2 or pi rad/samples
		vector<int32_t> filterCoeff{
			static_cast<int32_t>(0.009821820708929 * INT16_MAX),
			static_cast<int32_t>(-0.04103122121551 * INT16_MAX),
			static_cast<int32_t>(-0.05751222300691 * INT16_MAX),
			static_cast<int32_t>(0.01999189052168 * INT16_MAX),
			static_cast<int32_t>(0.1968631695768 * INT16_MAX),
			static_cast<int32_t>(0.3525450588422 * INT16_MAX),
			static_cast<int32_t>(0.3525450588422 * INT16_MAX),
			static_cast<int32_t>(0.1968631695768 * INT16_MAX),
			static_cast<int32_t>(0.01999189052168 * INT16_MAX),
			static_cast<int32_t>(-0.05751222300691 * INT16_MAX),
			static_cast<int32_t>(-0.04103122121551 * INT16_MAX),
			static_cast<int32_t>(0.009821820708929 * INT16_MAX)
		};


		FilterUpsamplingFir<complex<int32_t>, complex<int32_t>, complex<int32_t>, int32_t, ratio> filter(filterCoeff);

		// ---- Display of filter parameters
		cout << "Upsampling Ratio: " << ratio << '\n';
		for (size_t index = 0; index< filterCoeff.size(); ++index)
			cout << filterCoeff[index] << " ";
		cout << '\n';

		// ----- Impulse response
		cout << "Filter inpulse response" << '\n';
		vector<complex<int32_t> > input{ { 1, 1 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };
		vector<complex<int32_t> > output(input.size() * ratio);
		filter.reset();
		filter.step(input, output);
		cout << "Input:\n";
		for (size_t index = 0; index< input.size(); ++index)
			cout << input[index] << " ";
		cout << '\n';
		cout << "Output:\n";
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
		cout << '\n';


		// ----- Sinewave response
		cout << "Sinewave response" << '\n';
		GenSine<complex<int32_t> > gen(2.06, 10000);
		vector<complex<int32_t> > sineSignal(1000);
		gen.step(sineSignal);
		output.resize(sineSignal.size() * ratio);
		filter.reset();
		filter.step(sineSignal, output);
		string sineFilename{ "sine_input.txt" };
		string filterOutFilename{ "filter_output.txt" };
		ofstream osIn(sineFilename);
		ofstream osOut(filterOutFilename);
		saveAsciiSamples(output, osOut);
		saveAsciiSamples(sineSignal, osIn);
		cout << "Filter input saved in :" << sineFilename << '\n';
		cout << "Filter output saved in : " << filterOutFilename << '\n';

	}
	return false;
}

#endif

