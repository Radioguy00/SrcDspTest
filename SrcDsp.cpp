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
#include "coefficients.h"


#ifdef _WIN32
#include "..\..\SrcLibraries\SrcDsp\filters.h"
#include "..\..\SrcLibraries\SrcDsp\generators.h"
#include "..\..\SrcLibraries\SrcDsp\files.h"
#include "..\..\SrcLibraries\SrcDsp\upsampling_filters.h"
#include "..\..\SrcLibraries\SrcDsp\dnsampling_filters.h"
#include "..\..\SrcLibraries\SrcDsp\modulators.h"
#include "..\..\SrcLibraries\SrcDsp\correlators.h"
#include "..\..\SrcLibraries\SrcDsp\mixers.h"
#else
#include "../../src_libraries/dsp/filters.h"
#include "../../src_libraries/dsp/generators.h"
#include "../../src_libraries/dsp/files.h"
#include "../../src_libraries/dsp/upsampling_filters.h"
#include "../../src_libraries/dsp/dnsampling_filters.h"
#include "../../src_libraries/dsp/modulators.h"
#include "../../src_libraries/dsp/correlators.h"
#include "../../src_libraries/dsp/mixers.h"
#endif

/*-----------------------------------------------------------------------------
Test ModulatorSdpsk<T>
@tparam T Type of the output of the modulator

@param bitsFile File in which the bits used to modulate the waveform are stored
@param outFile File in which the output of the modulator is stored


------------------------------------------------------------------------------*/
template<class T>
bool testModulatorSdpsk(std::string bitsFile, std::string outFile)
{
	using namespace std;
	using namespace dsptl;
	const size_t nbrBits = 100; // Number of bits to run through the modulator

	// Create the bit pattern
	vector<int8_t> bits(nbrBits);
	// Create the modulator object
	ModulatorSdpsk<T> mod{};
	// Create the ouput of the modulator
	vector<complex<T> > out(bits.size());
	// Run the modulator
	mod.step(bits, out);
	// Save the output file
	ofstream osBits(bitsFile);
	ofstream osOut(outFile);
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
		mixer.setFrequency(0.05);
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

	GenSine<InType> gen(0.1, 31000);
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

//	testFixedPatternCorrelator<int16_t>();

//	testFiles();
//	testGenerators();
//	testFilters<complex<int16_t>, complex<int16_t>, complex<int32_t>, int32_t>();
	testMixers<std::complex<int16_t>, std::complex<int16_t>, int16_t, 4096 >();
	//	testUpsamplingFilters();
//	testModulatorSdpsk<float>("Input.txt", "Output.txt");
//	testModulatorSdpsk<double>("Input.txt", "Output.txt");
//	testModulatorSdpsk<int16_t>("Input.txt", "Output.txt");
//	testModulatorSdpsk<int32_t>("Input.txt", "Output.txt");

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


