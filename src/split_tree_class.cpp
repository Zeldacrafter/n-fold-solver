#include "template.h"

template <typename U>
class Node {
public:
    const static size_t NO_PARENT = std::numeric_limits<size_t>::max();

    U val;
    int childrenCnt;
    size_t parentIdx;

    Node(U _val, size_t pIdx) : val{_val}, parentIdx{pIdx}, childrenCnt{U(0)} {}
};

template <typename U>
class SplitTree {
public:
    std::vector<Node<U>> tree;
    std::vector<size_t> freeIndices;

    int add(U value, size_t parInd) {
        assert(SZ(tree) > parInd || parInd == Node<U>::NO_PARENT);

        if(parInd != Node<U>::NO_PARENT) {
            tree[parInd].childrenCnt++;
        }

        if(SZ(freeIndices)) {
            int idx = freeIndices.back();
            freeIndices.pop_back();
            tree[idx] = Node<U>(value, parInd);
            return idx;
        } else {
            tree.emplace_back(value, parInd);
            return SZ(tree) - 1;
        }
    }

    void remove(int index) {
        assert(SZ(tree) > index);
        freeIndices.push_back(index);
        if(tree[index].parentIdx != Node<U>::NO_PARENT) {
            int parIdx = tree[index].parentIdx;
            if(--tree[parIdx].childrenCnt == 0) {
                remove(parIdx);
            }
        }
    }

    std::vector<U> constructPath(int idx) {
        assert(idx < SZ(tree));

        std::vector<U> res;
        while(tree[idx].parentIdx != Node<U>::NO_PARENT) {
            res.push_back(tree[idx].val);
            idx = tree[idx].parentIdx;
        }
        std::reverse(ALL(res));

        return res;
    }

    void clear() {
        tree.clear();
        freeIndices.clear();
    }
};
