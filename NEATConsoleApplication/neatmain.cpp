#include "Population.h"
#include "experiments.h"

using namespace std;

int main(int argc, char *argv[]) {
	srand((unsigned)time(NULL));
	if (argc != 2) {
		cerr << "A NEAT parameters file (.ne file) is required to run the experiments!" << endl;

		getchar();
		return -1;
	}

	NEAT::loadNeatParams(argv[1], true);
	cout << "loaded" << endl;
	cout << "Press any key: ";
	getchar();

	NEAT::Population *p = 0;
	p = pole2TestRealTime();
	if (p)
		delete p;
	getchar();
	return 0;
}

