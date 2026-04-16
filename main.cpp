#include <iostream>
#include "Sequence.h"

using namespace std;

int main() {

    Sequence seq("ATGCGT");
    // Sequence seq("AAAAAA");
    // Sequence seq("ATATAT");
    // Sequence seq("GCGCGC");

    cout << "Original: " << seq.toString() << endl;

    // INSERT
    seq.insert(2, 'A');
    cout << "After insert: " << seq.toString() << endl;

    // DELETE
    seq.remove(3);
    cout << "After delete: " << seq.toString() << endl;

    // UPDATE
    seq.update(1, 'G');
    cout << "After update: " << seq.toString() << endl;

    // LENGTH
    cout << "Length: " << seq.length() << endl;

    return 0;
}