#ifndef GLOBALFUNCTION_HPP
#define GLOBALFUNCTION_HPP
#include "octreeNode.hpp"

octreeNode* buildInitOctree()
{
	octreeNode* root = new octreeNode();
	return root;
}
#endif
