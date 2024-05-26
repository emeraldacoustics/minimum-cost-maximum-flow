#include <queue>
#include <vector>

using namespace std;

const int inf = 0x7f7f7f7f;
const int maxsz = 1000;

class MinCostMaxFlowEdge
{
private:
protected:
public:
	int u;
	int v;
	int c;
	int cost;
	int f;

	MinCostMaxFlowEdge(void);
	MinCostMaxFlowEdge(const int &, const int &, const int &, const int &, const int &);
};

inline MinCostMaxFlowEdge::MinCostMaxFlowEdge(void) : u(0), v(0), c(0), cost(0), f(0)
{

}

inline MinCostMaxFlowEdge::MinCostMaxFlowEdge(const int & u, const int & v, const int & c, const int & cost, const int & f) : u(u), v(v), c(c), cost(cost), f(f)
{

}

int n;
int sz;
vector<MinCostMaxFlowEdge> edges;
vector<int> Eid[maxsz];
int a[maxsz];
int d[maxsz];
bool inq[maxsz];
int p[maxsz];

inline void initialize(void)
{
	int i;

	sz = n;
	edges.clear();
	for(i=0 ; i<sz ; ++i)
		Eid[i].clear();
}

inline void add_edge(const int & u, const int & v, const int & c, const int & cost)
{
	edges.push_back(MinCostMaxFlowEdge(u, v, c, cost, 0));
	Eid[u].push_back(edges.size()-1);
	edges.push_back(MinCostMaxFlowEdge(v, u, 0, -cost, 0));
	Eid[v].push_back(edges.size()-1);
}

inline bool bellman_ford(const int & s, const int & t, int & flow, int & cost)
{
	int i;
	int u;
	MinCostMaxFlowEdge * e;
	queue<int> que;

	memset(d, inf, sizeof(d[0])*sz);
	memset(inq, false, sizeof(inq[0])*sz);
	a[s] = inf;
	d[s] = 0;
	inq[s] = true;
	for(que.push(s) ; !que.empty() ; que.pop())
	{
		u = que.front();
		inq[u] = false;

		for(i=0 ; i<Eid[u].size() ; ++i)
		{
			e = &edges[Eid[u][i]];

			if(e->c>e->f && d[e->v]>d[u]+e->cost)
			{
				a[e->v] = min(a[u], e->c-e->f);
				d[e->v] = d[u]+e->cost;
				p[e->v] = Eid[u][i];
				if(!inq[e->v])
				{
					que.push(e->v);
					inq[e->v] = true;
				}
			}
		}
	}

	if(d[t] == inf)
		return false;
	else
	{
		flow += a[t];
		cost += a[t]*d[t];
		for(u=t ; u!=s ; u=edges[p[u]].u)
		{
			edges[p[u]].f += a[t];
			edges[p[u]^1].f -= a[t];
		}

		return true;
	}
}

inline int minimum_cost_maximum_flow(const int & s, const int & t, int & cost)
{
	int flow;

	flow = 0;
	cost = 0;
	for(; bellman_ford(s, t, flow, cost); );

	return flow;
}
