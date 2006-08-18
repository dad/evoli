#ifndef SEQUENCE_HH
#define SEQUENCE_HH

#include <vector>
#include "tools.hh"

using namespace std;

class Sequence {
private:
	Sequence();
protected:
	vector<int> m_sequence;
	bool m_modified;

	virtual ~Sequence(void) {}

public:
	Sequence(const vector<int>& v) { m_sequence = v; }
	Sequence(const int length) : m_sequence(length) {}

	virtual uint16 length() const { return m_sequence.size(); }
	bool modified() const { return m_modified; }
	int& operator[](const int index) {
		m_modified = true;
		return m_sequence[index];
	}

	int operator[](const int index) const { return m_sequence[index]; }
	void clear() {
		m_sequence.clear();
		m_modified = true;
	}

	typedef vector<int>::iterator iterator;
	typedef vector<int>::const_iterator const_iterator;

	iterator begin() { return m_sequence.begin(); }
	iterator end() { return m_sequence.end(); }

	const_iterator begin() const { return m_sequence.begin(); }
	const_iterator end() const { return m_sequence.end(); }

	bool operator==(const Sequence& s) const {
		bool identical = (length() == s.length());
		const_iterator qit = begin();
		const_iterator pit = s.begin();
		for (; qit != end() && pit != s.end() && identical; qit++, pit++) {
			identical = (*pit == *qit);
		}
		return identical;
	}
};

#endif //SEQUENCE_HH
