#include <queue>
#include <vector>

using namespace std;

const int inf = 0x7f7f7f7f;
const long long infll = 0x7f7f7f7f7f7f7f7fll;
const long double infl = 1e20;
const long double eps = 1e-12;

template <class T>
inline int sgn(const T & x)
{
	return (x > eps) - (x < -eps);
}

template <class T, class U = T>
class MinimumCostMaximumFlow
{
public:
	static const int maxn = maxN * 2 + 2;

	class Edge
	{
	public:
		int u, v;
		T c, f;
		U cost;

		Edge(void) : u(0), v(0), c(0), cost(0), f(0)
		{

		}

		Edge(const int & u, const int & v, const T & c, const U & cost, const T & f) : u(u), v(v), c(c), cost(cost), f(f)
		{

		}
	};

	int n;
	vector<Edge> edges;
	vector<int> Eid[maxn];
	T a[maxn];
	U d[maxn];
	bool inq[maxn];
	int p[maxn];

	MinimumCostMaximumFlow(void) : n(0)
	{

	}

	void resize(const int & n)
	{
		edges.clear();
		for (int i = 0; i < this->n; i++)
			Eid[i].clear();
		this->n = n;
	}

	void add_edge(const int & u, const int & v, const T & c, const U & cost)
	{
		edges.push_back(Edge(u, v, c, cost, 0));
		Eid[u].push_back(edges.size() - 1);
		edges.push_back(Edge(v, u, 0, -cost, 0));
		Eid[v].push_back(edges.size() - 1);
	}

	bool bellman_ford(const int & s, const int & t, T & flow, U & cost)
	{
		queue<int> que;

		/**/
		a[s] = inf;
		fill(d, d + n, infl);
		d[s] = 0;
		memset(inq, false, sizeof inq[0] * n);
		inq[s] = true;
		for (que.push(s); !que.empty(); que.pop())
		{
			const int & u = que.front();
			inq[u] = false;

			for (const auto & id : Eid[u])
			{
				Edge & e = edges[id];
				if (sgn(e.c - e.f) > 0 && sgn(d[e.v] - (d[u] + e.cost)) > 0)
				{
					a[e.v] = min(a[u], e.c - e.f);
					d[e.v] = d[u] + e.cost;
					p[e.v] = id;
					if (!inq[e.v])
					{
						inq[e.v] = true;
						que.push(e.v);
					}
				}
			}
		}

		/**/
		if (sgn(d[t] - infl) == 0)
			return false;
		else
		{
			flow += a[t];
			cost += d[t] * a[t];
			for (int u = t; u != s; u = edges[p[u]].u)
			{
				edges[p[u]].f += a[t];
				edges[p[u] ^ 1].f -= a[t];
			}
			return true;
		}
	}

	T minimum_cost_maximum_flow(const int & s, const int & t, U & cost)
	{
		T flow = 0;
		cost = 0;
		for (; bellman_ford(s, t, flow, cost); );
		return flow;
	}
};
