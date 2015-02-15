// SrcDsp.cpp : Defines the entry point for the console application.
//

#ifdef _WIN32
#include "stdafx.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

#ifdef _WIN32
#include "..\..\SrcLibraries\SrcDsp\filters.h"
#include "..\..\SrcLibraries\SrcDsp\generators.h"
#include "..\..\SrcLibraries\SrcDsp\files.h"
#include "..\..\SrcLibraries\SrcDsp\upsampling_filters.h"
#include "..\..\SrcLibraries\SrcDsp\modulators.h"
#else
#include "../../src_libraries/dsp/filters.h"
#include "../../src_libraries/dsp/generators.h"
#include "../../src_libraries/dsp/files.h"
#include "../../src_libraries/dsp/upsampling_filters.h"
#include "../../src_libraries/dsp/modulators.h"
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



bool testFilters();
bool testGenerators();
bool testFiles();
bool testUpsamplingFilters();
bool testModulators();
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
	testFiles();
	testGenerators();
	testFilters();
	testUpsamplingFilters();
	testModulatorSdpsk<float>("Input.txt", "Output.txt");
	testModulatorSdpsk<double>("Input.txt", "Output.txt");
	testModulatorSdpsk<int16_t>("Input.txt", "Output.txt");
	testModulatorSdpsk<int32_t>("Input.txt", "Output.txt");

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

bool testFilters()
{
	using namespace std;

	cout << "Program is starting" << '\n';
	cout << "Floating point test" << '\n';

	{
		cout << "Creation of filter object" << '\n';
		vector<float> filterCoeff{ 1.34f, 4.56f, 0.34f };
		FilterFir<float, double, double, float> filter(filterCoeff);

		cout << "Filter inpulse response" << '\n';
		vector<float> input{ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		vector<double> output(input.size());
		filter.step(input, output);
#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif
		cout << '\n';

		cout << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif
		cout << '\n';
	}

	cout << '\n' << "Integer test" << '\n';
	{
		cout << "Creation of filter object" << '\n';
		vector<int16_t> filterCoeff{ 234, -23, 12 };
		FilterFir<int16_t, int32_t, int32_t, int16_t> filter(filterCoeff);

		cout << "Filter inpulse response" << '\n';
		vector<int16_t> input{ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		vector<int32_t> output(input.size());
		filter.step(input, output);

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif

		cout << '\n' << "Delayed Filter inpulse response" << '\n';
		input = { 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 };
		output.resize(input.size());
		filter.step(input, output);

#ifdef CPLUSPLUS11
		for (const auto& elt : output)
			cout << elt << " ";
#else
		for (size_t index = 0; index< output.size(); ++index)
			cout << output[index] << " ";
#endif
		cout << '\n';
	}

	return false;
}


bool testModulators()
{



	return false;
}
