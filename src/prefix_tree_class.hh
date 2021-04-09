#ifndef N_FOLD_PREFIX_TREE_CLASS_HH
#define N_FOLD_PREFIX_TREE_CLASS_HH

#include <iostream>

#include "utils.hh"

/**
 * Data structure for saving multiple possibly overlapping sequences efficiently.
 * @tparam U The type of the values in the sequences.
 */
template <typename U>
class prefix_tree {
  private:
    struct Node {
      public:
        U val;
        size_t parentIdx;
        int childrenCnt;

        Node(U _val, size_t pIdx) : val{_val}, parentIdx{pIdx}, childrenCnt{U(0)} {}
    };

    std::vector<Node> tree;
    std::vector<size_t> freeIndices;

  public:
    const static size_t NO_PARENT = std::numeric_limits<size_t>::max();

    /**
     * Add an element to the tree structure.
     * @param value The value of the element to add.
     * @param parInd The index of the parent of the added element.
     * @return The index of the addes element.
     */
    int add(U value, size_t parInd) {
        assert(tree.size() > parInd || parInd == NO_PARENT);

        if(parInd != NO_PARENT) {
            tree[parInd].childrenCnt++;
        }

        if(freeIndices.size()) {
            int idx = freeIndices.back();
            freeIndices.pop_back();
            tree[idx] = Node(value, parInd);
            return idx;
        } else {
            tree.emplace_back(value, parInd);
            return tree.size() - 1;
        }
    }

    /**
     * Remove an element from the tree structure.
     * @param index The index of the element to remove.
     * @param removeParents Flag to indicate whether parents should be recursively removed if they have no children anymore.
     *                      true by default but should be set to false if new children might be added to the parent after this
     *                      remove operation.
     */
    void remove(size_t index, bool removeParents = true) {
        assert(tree.size() > index);
        freeIndices.push_back(index);
        if(tree[index].parentIdx != NO_PARENT) {
            int parIdx = tree[index].parentIdx;
            if(--tree[parIdx].childrenCnt == 0 && removeParents) {
                remove(parIdx);
            }
        }
    }

    /**
     * Construct the path with all values from the root to the specified element.
     * @param index The index of the last node in the path.
     * @return Vector of values on that path.
     */
    std::vector<U> constructPath(size_t index) {
        assert(index < tree.size());

        std::vector<U> res;
        while(tree[index].parentIdx != NO_PARENT) {
            res.push_back(tree[index].val);
            index = tree[index].parentIdx;
        }
        std::reverse(res.begin(), res.end());

        return res;
    }

    /**
     * Remove all elements from the tree. This keeps memory allocated so that it can be reused.
     */
    void clear() {
        tree.clear();
        freeIndices.clear();
    }
};

#endif //N_FOLD_PREFIX_TREE_CLASS_HH