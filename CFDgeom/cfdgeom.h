/*
CFDgeom
	1|	CFD-ready mesh info archiving
	2|	Cell interconnectivity and labeling
*/

#ifndef CFDGEOM_H
#define CFDGEOM_H

#include<iostream>
#include<vector>
#include<map>
#include<string>
#include<memory>
#include<utility> // std::pair
#include<Eigen/Dense>
#include<mshio/mshio.h>
#include"utilgeom.h"

typedef unsigned int uint;
typedef Eigen::Vector3d coor;

namespace cfdgeom
{
struct fcontain // face container
{
	std::vector<uint> fnode;
	std::string fbound;
};

struct finfo // face info container
{
	coor fcentroid;
	coor fnormal; // e
	double farea;
	std::string bound;
};

struct ccontain // cell container
{
	std::vector<uint> cface;
	std::string cdomain;
};

struct cinfo // cell info container
{
	coor ccentroid;
	double cvolume;
	std::string domain;
};

class mesh
{	
	private:
	void gen_container(std::map<uint, coor>*, std::map<uint, fcontain>*, std::map<uint, ccontain>*, mshio::MshSpec);
	void gen_neighbor(std::map<uint, uint[2]>*, std::map<uint, uint>*, std::map<uint, fcontain>*, std::map<uint, ccontain>*);
	void gen_info  (std::map<uint, finfo>*, std::map<uint, cinfo>*, std::map<uint, coor>*,
					std::map<uint, fcontain>*, std::map<uint, ccontain>*);

	public:
	std::unique_ptr<std::map<uint, finfo>> face;
	std::unique_ptr<std::map<uint, cinfo>> cell;
	std::unique_ptr<std::map<uint, uint[2]>> ccneigh; // face id, < cell 1 id, cell 2 id>
	std::unique_ptr<std::map<uint, uint>> fcneigh; // face id, cell id

	mesh(mshio::MshSpec spec)
	{
		std::map<uint, coor>* lnode = new std::map<uint, coor>;
		std::map<uint, fcontain>* lface = new std::map<uint, fcontain>;
		std::map<uint, ccontain>* lcell = new std::map<uint, ccontain>;
		
		// fill empty pointers
		face.reset(new std::map<uint, finfo>);
		cell.reset(new std::map<uint, cinfo>);
		ccneigh.reset(new std::map<uint, uint[2]>);
		fcneigh.reset(new std::map<uint, uint>)

		gen_container(lnode, lface, lcell, spec);
		gen_neighbor(ccneigh, fcneigh, lface, lcell);
		gen_info(face, cell, lnode, lface, lcell);

		// garbage collector
		delete lnode;
		delete lface;
		delete lcell;

		cout << "mesh obj. generated..." << endl;
	};

}; // mesh

void mesh::gen_container(plnode, plface, plcell, spec)
{
	/* generate mshio containers */
	
	// map out boundary (surface and volume) names by tag
	std::map<uint, std::string> tag_surf;
	std::map<uint, std::string> tag_vol;

	for(const auto& group: spec.physical_group)
	{
		switch(group.dim)
		{
		case 2:
		tag_surf.insert({group.tag, group.name});

		case 3:
		tag_vol.insert({group.tag, group.name});

		default:
		cout << "Invalid dimemsion size. Abort..." << endl
		return
		};
	}; // for

	// map out all node coor
	for(std::size_t i = 0; i < spec.nodes.num_entity_blocks; i++)
	{
		mshio::NodeBlock& block = spec.nodes.num_entity_blocks[i];

		for(std::size_t j = 0; j < block.num_nodes_in_block; j++)
		{
			coor node;

			for(size_t k = 0; k < 3; k++)
			{
				node(k) = (block.data[j * (3) + k]);
			};

			*plnode.insert({block.tags[j], node});
		};
	};

	// face and cell iter
	for(std::size_t i = 0; i < spec.elements.num_entity_blocks; i++)
	{
		mshio::ElementBlock& block = spec.elements.entity_blocks[i];

		switch(block.entity_dim)
		{
			// 2d always comes before 3d
			case 2:
			uint tag = block.entity_tag;
			uint type = block.element_tag;

			for(std::size_t j = 0; j < block.num_elements_in_block; j++)
			{
				fcontain newFace;
				newFace.bound = check_tag(tag_surf, tag); // utilgeom.h

				// append node tag list
				std::size_t n = mshio::nodes_per_element(type);
				std::vector<uint> nodes;

				for (std::size_t k = 0; k <= n; k++)
				{
					nodes.push_back(block.data[j * (n + 1) + k]);
				};

				newFace.fnode = nodes;
				*plface.insert({block.data[j * (n + 1)], newFace});

			}; // for

			case 3:
			uint tag = block.entity_tag;
			uint type = block.element_type;

			for(size_t j = 0; j < block.num_elements_in_block; j++)
			{
				ccontain newCell;
				newCell.domain = check_tag(tag_vol, tag);
				std::size_t n = mshio::nodes_per_element(type);
				std::vector<uint> nodes;

				for(std::size_t k = 0; k <= n; k++)
				{
					nodes.push_back(block.data[j * (n + 1) + k]);
				};

				// cell nodes to face ids
				std::vector<uint> lcface;

				for(std::pair<K, V> entry: *pface)
				{
					{
					if(match_vec(entry.second->fnode, nodes))
					{
						lcface.push_back(entry.first);
					};
					else if(lcface.size() == nodes.size())
					{
						break;
					};
					} // if else

				};

				newCell.cface = lcface;
				*plcell.insert({block.data[j * (n + 1)], newCell});
			}; // for

			default:
			cout << "Invalid dimemsion size. Abort..." << endl
			return
		}; // switch
	}; // for

}; // gen_container

void mesh::gen_neighbor(pccneigh, pfcneigh, plface, plcell)
{
	/* generate neighbors */
	// structured, a face can only be connected once (source) or twice (neighbor)

	for(std::pair<K, V> entry1: *plface)
	{
		switch(entry1.second->fbound)
		{
			case "none": // ccneigh
			uint newNeigh [2];
			for(std::pair<K, V> entry2: *plcell)
			{
				unsigned short int ctrl = 0;
				std::vector<uint> ref{entry1.first};
				{
				if(ctrl != 2 && match_vec(ref, entry2.second->cface))
				{
					newNeigh[ctrl] = entry2->first;
					ctrl++;
				};
				else if(ctrl == 2)
				{
					*pccneigh.insert({entry1.first, newNeigh});
					break;
				};
				}; // if else

			}; // "none"
		
			default: // fcneigh
			for(std::pair<K, V> entry2: *plcell)
			{
				std::vector<uint> ref{entry.first};
				
				if(match_vec(ref, entry2.second->cface))
				{
					*pfcneigh.insert({entry1.first, entry2.first});
					break;
				};

			}; // default

		}; // switch
	};	// for

}; // gen_neighbor

void mesh::gen_info(pface, pcell, plnode, plface, plcell)
{
	/* generate info */
	// face
	for(std::pair<K, V> entry: *plface)
	{
		finfo newFace;
		newFace.bound = entry.second->fbound;
		std::vector<coor> fcoor;

		for(auto i = entry.second->fcoor.begin(); i != entry.second->fcoor.end(); i++)
		{
			fcoor.push_back(*plnode[i]);
		};

		std::pair<std::pair<coor, coor>, double> info = info_2d(fcoor);

		newFace.fcentroid = info.first.first;
		newFace.fnormal = info.first.second;
		newFace.farea = info.second;

		*pface.insert({entry.first, newFace});
	}; // for

	// cell
	for(std::pair<K, V> entry: *plcell)
	{
		cinfo newCell;
		newCell.domain = entry->second.cdomain;

		std::pair<coor, double> info = info_3d(entry.second->cface, pface);

		newCell.ccentroid = info_3d.first;
		newCell.cvolume = info_3d.second;
		
		*pcell.insert({entry.first, newCell});
	}; // for

}; // gen_info

}; // namespace cfdgeom

#endif