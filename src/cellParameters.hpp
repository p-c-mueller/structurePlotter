s/*
 * cellParameters.hpp
 *
 *  Created on: Jul 11, 2023
 *      Author: pemueller
 */

#ifndef SRC_CELLPARAMETERS_HPP_
#define SRC_CELLPARAMETERS_HPP_

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <utility>

using namespace Eigen;
using namespace std;

class Atom {
public:
	Atom(){}
	Atom(string element, Vector3d coordinatesFrac, Vector3d coordinatesReal, Vector3d rgb = Vector3d::Zero())
	{
		_element = element;
		_coordinatesFrac = coordinatesFrac;
		_coordinatesReal = coordinatesReal;
		_translation(0) = (int)coordinatesFrac(0);
		_translation(1) = (int)coordinatesFrac(1);
		_translation(2) = (int)coordinatesFrac(2);

		_rgb = rgb;
		_size = 0;
		_index = -1;
	}
	~Atom(){}

	void setRGB( Vector3d rgb ){ _rgb = rgb; }
	void setSize( double size ){ _size = size; }
	void setCoordinatesFrac ( Vector3d coordinates) { _coordinatesFrac = coordinates; }
	void setCoordinatesReal ( Vector3d coordinates) { _coordinatesReal = coordinates; }

	string getElement() { return _element; }
	Vector3d getCoordinatesFrac(){ return _coordinatesFrac; }
	Vector3d getCoordinatesReal(){ return _coordinatesReal; }
	Vector3d getTranslation(){ return _translation; }
	Vector3d getRGB(){ return _rgb; }
	double getSize(){ return _size; }

	void setIndex( int index ){ _index = index; }
	int getIndex(){ return _index; }

	bool operator == (const Atom other) const
	{
		return (_element == other._element && (_coordinatesFrac - other._coordinatesFrac).norm() < 0.0001);
	}
private:
	string _element;
	Vector3d _coordinatesFrac;
	Vector3d _coordinatesReal;
	Vector3d _translation;
	Vector3d _rgb;
	double _size;
	int _index;
};
typedef boost::shared_ptr<Atom> AtomPtr;

class Bond {
public:
	Bond(){}

	Bond( pair<AtomPtr,AtomPtr> atoms, double size, Vector3d rgb = Vector3d::Zero())
	{
		_atoms = atoms;
		_length = (atoms.first->getCoordinatesReal() - atoms.second->getCoordinatesReal()).norm();
		_size = size;
		_rgb = rgb;
	}
	~Bond(){}

	AtomPtr getFrom(){ return _atoms.first; }
	AtomPtr getTo(){ return _atoms.second; }
	double getLength(){ return _length; }
	double getSize() {return _size;}
	Vector3d getRGB(){ return _rgb; }

	bool operator == (const Bond other) const
	{
		return ( ((*_atoms.first == *other._atoms.first && *_atoms.second == *other._atoms.second)
				|| (*_atoms.first == *other._atoms.second && *_atoms.second == *other._atoms.first)));
	}

	Vector3d getCoordinatesFrac()
	{
		if ( _atoms.first->getCoordinatesFrac()(2) < _atoms.second->getCoordinatesFrac()(2) )
			return _atoms.first->getCoordinatesFrac();
		else
			return _atoms.second->getCoordinatesFrac();
	}

	Vector3d getCoordinatesReal()
	{
		if ( _atoms.first->getCoordinatesReal()(2) >= _atoms.second->getCoordinatesReal()(2) )
			return _atoms.first->getCoordinatesReal();
		else
			return _atoms.second->getCoordinatesReal();
	}

	Vector3d getEnd()
	{
		if ( _atoms.first->getCoordinatesReal()(2) <= _atoms.second->getCoordinatesReal()(2) )
			return _atoms.first->getCoordinatesReal();
		else
			return _atoms.second->getCoordinatesReal();
	}

	Vector3d getDirection()
	{
		//if ( _atoms.first->getCoordinatesReal()(2) < _atoms.second->getCoordinatesReal()(2) )
			//return _atoms.first->getCoordinatesReal() - _atoms.second->getCoordinatesReal();
		//else
			return _atoms.second->getCoordinatesReal() - _atoms.first->getCoordinatesReal();
	}

private:
	pair<AtomPtr,AtomPtr> _atoms;
	double _length;
	double _size;
	Vector3d _rgb;
};
typedef boost::shared_ptr<Bond> BondPtr;

class Polyhedron
{};
typedef boost::shared_ptr<Polyhedron> PolyPtr;

class BondPreset{
public:
	BondPreset(){}
	BondPreset( pair<string,string> types, pair<double,double> length, Vector3d rgb, double size, int boundary )
	{
		_types = types;
		_length = length;
		_rgb = rgb;
		_size = size;
		_boundary = boundary;
	}

	virtual ~BondPreset(){}

	pair<string,string> getTypes(){return _types;}
	pair<double,double> getLength(){ return _length;}
	Vector3d getRGB(){return _rgb;}
	double getSize(){ return _size;}
	int getBoundary(){ return _boundary;}

private:
//	UnitCellPtr _unitCell;
	pair<string,string> _types;
	pair<double,double> _length;
	Vector3d _rgb;
	double _size;
	int _boundary;
};
typedef boost::shared_ptr<BondPreset> BondPresetPtr;

class UnitCell{
public:
	UnitCell(){}
	~UnitCell(){}

	void setLattice( Matrix3d lattice ) { _lattice = lattice; }
	MatrixXd getLattice() { return _lattice; }
	void setBoundary( pair<Vector3d,Vector3d> boundary ) { _boundary = boundary; }
	pair<Vector3d,Vector3d> getBoundary(){ return _boundary; }

	void addAtom( AtomPtr atom )
	{
		int index = _atoms.size();
		_atoms.push_back(atom);
		_atoms[index]->setIndex(index);
	}

	AtomPtr getAtom( AtomPtr atom )
	{
		for ( auto it : _atoms )
		{
			if ( *atom == *it )
				return it;
		}
	}

	AtomPtr getAtomAtPosition( Vector3d position )
	{
		int index = -1;
		double distMin = 0;

		for ( int i = 0; i < _atoms.size(); ++i )
		{
			double dist = (_atoms[i]->getCoordinatesReal() - position).norm();
			if( (index < 0 || dist < distMin) && _atoms[i]->getElement() != "NULL" )
			{
				index = i;
				distMin = dist;
			}
		}

		return _atoms[index];
	}

	bool atomExists( AtomPtr atom )
	{
		for ( auto it : _atoms )
			if ( *it == *atom )
				return true;
		return false;
	}

	vector<AtomPtr> getAtoms(){ return _atoms; }
	void setOrientation( Matrix3d orientation ){ _orientation = orientation; }
	Matrix3d getOrientation(){ return _orientation; }

	void addBondPreset( BondPresetPtr preset ){ _bondPresets.push_back(preset); }
	vector<BondPresetPtr> getBondPresets(){ return _bondPresets; }

	void addBond( BondPtr bond )
	{
		if ( bond->getLength() == 0 )
			return;

		for ( auto _bond : _bonds )
		{
			if ( *bond == *_bond )
			{
				return;
			}
		}

		_bonds.push_back(bond);
	}
	vector<BondPtr> getBonds(){return _bonds;}

private:

	Matrix3d _lattice;
	pair<Vector3d,Vector3d> _boundary;
	Matrix3d _orientation;
	vector<AtomPtr> _atoms;
	vector<BondPtr> _bonds;
	// vector of bond presets
	vector<BondPresetPtr> _bondPresets;
};
typedef boost::shared_ptr<UnitCell> UnitCellPtr;


class Object
{
public:
	Object(){};
	~Object(){};
	Object( string type, Vector2d position, Vector3d rgb, double z, double size = 0, Vector2d direction = Vector2d::Zero() )
	{
		_type = type;
		_position = position;
		if ( _type == "line" )
			_z = z;
		else
			_z = z + size;

		_size = size;
		_direction = direction;
		_rgb = rgb;
	}

	string getType(){return _type;}
	double getSize(){return _size;}
	Vector2d getPosition(){ return _position; }
	Vector2d getDirection(){ return _direction; }
	Vector3d getRGB(){return _rgb;}
	void addX( double x ) { _position(0) += x; _direction(0)+= x; }
	void addY( double y ) { _position(1) += y; _direction(1)+=y; }

	bool operator < ( const Object other ) const
		{return _z < other._z;}
	bool operator < ( const boost::shared_ptr<Object> other ) const
		{return _z < other->_z;}

	double getZ(){return _z;}


private:
	string _type;
	Vector2d _position;
	double _z;
	double _size;
	Vector2d _direction;
	Vector3d _rgb;

};
typedef boost::shared_ptr<Object> ObjectPtr;

#endif /* SRC_CELLPARAMETERS_HPP_ */
