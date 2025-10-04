#pragma once

struct edge {
    int cost;
    int from;
    int where_;
};

struct compare {
    bool operator() (edge const& a, edge const& b) const {
        return a.cost > b.cost;
    }
};

struct comparenum {
    bool operator() (edge const& a, edge const& b) const {
        return a.from < b.from;
    }
};
