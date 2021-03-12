#include "template.h"
#include "bruteforce.cpp"

int main() {
    using namespace std;

    int n, r, s, t;
    cin >> n >> r >> s >> t;

    NFold<int> nfold(n, r, s, t);
    cin >> nfold;

    cout << nfold << endl;

    cout << "A_init:" << endl
         << constructAInit(nfold) << endl;

    //cout << bruteForce(nfold) << endl;
    return 0;
}