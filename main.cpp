#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <array>
#include <string>
#include <utility>
using namespace std;

static const char EPS = '\0'; ///< 表示ε空转换的特殊字符

//-------------------- Regex Parser --------------------

/**
 * @brief 正则表达式的AST节点。
 *
 * RegexNode用于表示正则表达式的抽象语法树(AST)节点，其中包含了四种类型:
 * - CHAR: 单字符
 * - CONCAT: 连接运算
 * - UNION: 并运算(使用'+'表示)
 * - STAR: Kleene闭包运算('*')
 */
struct RegexNode {
    /**
     * @brief 表示RegexNode的类型枚举。
     */
    enum Type { CHAR, CONCAT, UNION, STAR } type; ///< 节点类型
    char ch; ///< 对于CHAR类型，表示此节点的字符；其它类型此字段无意义
    RegexNode *left;  ///< 左子节点
    RegexNode *right; ///< 右子节点

    /**
     * @brief 构造RegexNode节点。
     * @param t 节点类型
     * @param c 字符(仅CHAR类型使用)
     * @param l 左子节点指针
     * @param r 右子节点指针
     */
    RegexNode(Type t, char c = 0, RegexNode *l = nullptr, RegexNode *r = nullptr)
        : type(t), ch(c), left(l), right(r) {}
};

/**
 * @brief 正则表达式解析器(使用递归下降法)。
 *
 * 从输入的正则表达式字符串中构建抽象语法树(AST)。
 * 支持字符'0','1'，以及'+'表示并、'*'表示Kleene星闭包、'()'分组和隐式连接。
 */
class RegexParser {
public:
    /**
     * @brief 构造函数
     * @param s 正则表达式字符串
     */
    explicit RegexParser(string s) : str(std::move(s)), pos(0) {}

    /**
     * @brief 解析入口函数。
     * @return 构造完成的AST根节点指针
     */
    RegexNode *parse() {
        return parseUnion();
    }

private:
    string str; ///< 待解析的正则表达式字符串
    int pos;    ///< 当前解析到的位置

    /**
     * @brief 解析Union(并运算：'+')
     * @return 语法子树根节点
     */
    RegexNode *parseUnion() {
        RegexNode *node = parseConcat();
        while (pos < (int) str.size() && str[pos] == '+') {
            pos++;
            RegexNode *right = parseConcat();
            node = new RegexNode(RegexNode::UNION, 0, node, right);
        }
        return node;
    }

    /**
     * @brief 解析Concat(连接运算)
     * @return 语法子树根节点
     */
    RegexNode *parseConcat() {
        RegexNode *node = parseStar();
        while (pos < (int) str.size() && (str[pos] != ')' && str[pos] != '+')) {
            RegexNode *right = parseStar();
            node = new RegexNode(RegexNode::CONCAT, 0, node, right);
        }
        return node;
    }

    /**
     * @brief 解析Kleene星闭包
     * @return 语法子树根节点
     */
    RegexNode *parseStar() {
        RegexNode *node = parseBase();
        while (pos < (int) str.size() && str[pos] == '*') {
            pos++;
            node = new RegexNode(RegexNode::STAR, 0, node, nullptr);
        }
        return node;
    }

    /**
     * @brief 解析基本单元：字符或括号表达式
     * @return 语法子树根节点
     */
    RegexNode *parseBase() {
        if (pos >= (int) str.size()) return nullptr;
        if (str[pos] == '(') {
            pos++;
            RegexNode *node = parseUnion();
            if (pos < (int) str.size() && str[pos] == ')') {
                pos++;
            }
            return node;
        } else if (str[pos] == '0' || str[pos] == '1') {
            RegexNode *node = new RegexNode(RegexNode::CHAR, str[pos]);
            pos++;
            return node;
        }
        return nullptr;
    }
};

//-------------------- NFA定义 --------------------

/**
 * @brief 非确定性有限自动机(NFA)定义结构。
 *
 * 使用Thompson构造法将正则表达式转换为ε-NFA，再经过ε-消除等步骤。
 */
struct NFA {
    /**
     * @brief NFA状态结构。
     */
    struct State {
        int id; ///< 状态编号
        bool accept = false; ///< 是否是接受态
        map<char, vector<int> > trans; ///< 转移函数，键为输入字符，值为目标状态集合
    };

    vector<State> states; ///< 状态集合
    int start; ///< 起始状态ID

    /**
     * @brief 创建一个新的NFA状态
     * @param accept 此状态是否是接受态
     * @return 新状态的ID
     */
    int new_state(bool accept = false) {
        State s;
        s.id = (int) states.size();
        s.accept = accept;
        states.push_back(s);
        return s.id;
    }
};

/**
 * @brief 一个NFA片段，用于Thompson构造中间步骤。
 */
struct NFAFragment {
    int start;  ///< 片段的起始状态
    int accept; ///< 片段的接受状态(仅一个接受状态)
};

/**
 * @brief Thompson构造类，将正则的AST构造为ε-NFA。
 */
class Thompson {
public:
    /**
     * @brief 构造函数
     * @param root 正则表达式AST根节点
     */
    Thompson(RegexNode *root) : r(root) {}

    /**
     * @brief 使用Thompson构造法将AST转为ε-NFA
     * @return 构造好的NFA
     */
    NFA build() {
        nfa = NFA();
        NFAFragment frag = buildFragment(r);
        nfa.states[frag.accept].accept = true;
        nfa.start = frag.start;
        return nfa;
    }

private:
    RegexNode *r; ///< AST根节点
    NFA nfa;       ///< 构造中的NFA

    /**
     * @brief 添加一个状态转移
     * @param from 源状态
     * @param to 目标状态
     * @param c 转移字符(或EPS)
     */
    void add_transition(int from, int to, char c) {
        nfa.states[from].trans[c].push_back(to);
    }

    /**
     * @brief 针对特定的RegexNode构造NFA片段
     * @param node AST节点
     * @return 对应的NFAFragment
     */
    NFAFragment buildFragment(RegexNode *node) {
        switch (node->type) {
            case RegexNode::CHAR: {
                int s = nfa.new_state();
                int a = nfa.new_state();
                add_transition(s, a, node->ch);
                return {s, a};
            }
            case RegexNode::CONCAT: {
                NFAFragment f1 = buildFragment(node->left);
                NFAFragment f2 = buildFragment(node->right);
                add_transition(f1.accept, f2.start, EPS);
                return {f1.start, f2.accept};
            }
            case RegexNode::UNION: {
                int s = nfa.new_state();
                int a = nfa.new_state();
                NFAFragment f1 = buildFragment(node->left);
                NFAFragment f2 = buildFragment(node->right);
                add_transition(s, f1.start, EPS);
                add_transition(s, f2.start, EPS);
                add_transition(f1.accept, a, EPS);
                add_transition(f2.accept, a, EPS);
                return {s, a};
            }
            case RegexNode::STAR: {
                int s = nfa.new_state();
                int a = nfa.new_state();
                NFAFragment f = buildFragment(node->left);
                add_transition(s, f.start, EPS);
                add_transition(f.accept, a, EPS);
                add_transition(s, a, EPS);
                add_transition(f.accept, f.start, EPS);
                return {s, a};
            }
        }
        return {-1, -1};
    }
};

/**
 * @brief 消除NFA中的ε转移，得到无ε的NFA。
 */
class EpsilonRemover {
public:
    /**
     * @brief 构造函数
     * @param n 输入的ε-NFA
     */
    EpsilonRemover(const NFA &n) : infa(n) {}

    /**
     * @brief 消除ε转移
     * @return 无ε转移的NFA
     */
    NFA remove() {
        eps_closures.resize(infa.states.size());
        for (int i = 0; i < (int) infa.states.size(); i++) {
            compute_epsilon_closure(i);
        }

        NFA onfa;
        onfa.states.resize(infa.states.size());
        for (int i = 0; i < (int) infa.states.size(); i++) {
            onfa.states[i].id = i;
        }
        onfa.start = infa.start;

        for (int i = 0; i < (int) infa.states.size(); i++) {
            set<int> closure = eps_closures[i];
            bool is_accept = false;
            map<char, set<int> > combined;
            for (int cst: closure) {
                if (infa.states[cst].accept) is_accept = true;
                for (auto &kv: infa.states[cst].trans) {
                    char c = kv.first;
                    if (c == EPS) continue;
                    for (int nxt: kv.second) {
                        for (int ec: eps_closures[nxt]) {
                            combined[c].insert(ec);
                        }
                    }
                }
            }
            onfa.states[i].accept = is_accept;
            for (auto &kv: combined) {
                onfa.states[i].trans[kv.first] = vector<int>(kv.second.begin(), kv.second.end());
            }
        }

        return onfa;
    }

private:
    const NFA &infa; ///< 输入的ε-NFA
    vector<set<int> > eps_closures; ///< 每个状态的ε闭包集合

    /**
     * @brief 计算给定状态s的ε闭包
     * @param s 状态ID
     */
    void compute_epsilon_closure(int s) {
        if (!eps_closures[s].empty()) return;
        stack<int> st;
        st.push(s);
        eps_closures[s].insert(s);
        while (!st.empty()) {
            int u = st.top();
            st.pop();
            auto it = infa.states[u].trans.find(EPS);
            if (it != infa.states[u].trans.end()) {
                for (int nxt: it->second) {
                    if (eps_closures[s].find(nxt) == eps_closures[s].end()) {
                        eps_closures[s].insert(nxt);
                        st.push(nxt);
                    }
                }
            }
        }
    }
};

/**
 * @brief 确定性有限自动机(DFA)结构。
 */
struct DFA {
    /**
     * @brief DFA状态结构。
     */
    struct State {
        int id;       ///< 状态编号
        bool accept;  ///< 是否为接受态
        int t0, t1;   ///< 输入'0','1'时的转移
    };

    vector<State> states; ///< 状态集合
    int start; ///< 起始状态ID
    int trap;  ///< 陷阱态ID(无效转移的收容态)

    DFA() : start(-1), trap(-1) {}
};

/**
 * @brief 使用子集构造法(NFA->DFA)的类。
 */
class SubsetConstruction {
public:
    /**
     * @brief 构造函数
     * @param n 无ε的NFA
     */
    SubsetConstruction(const NFA &n) : infa(n) {}

    /**
     * @brief 将无ε NFA转换为DFA
     * @return 得到的DFA
     */
    DFA convert() {
        map<set<int>, int> state_map;
        vector<set<int> > dfa_sets;
        auto start_set = set<int>({infa.start});

        queue<set<int> > q;
        q.push(start_set);
        state_map[start_set] = 0;
        dfa_sets.push_back(start_set);

        struct TempState {
            int id;
            bool acc;
            int t0;
            int t1;
        };
        vector<TempState> tmp;

        while (!q.empty()) {
            auto cur = q.front();
            q.pop();
            int cid = state_map[cur];

            bool is_accept = false;
            for (auto s: cur) {
                if (infa.states[s].accept) {
                    is_accept = true;
                    break;
                }
            }

            auto go0 = move_set(cur, '0');
            auto go1 = move_set(cur, '1');

            int go0_id = get_state_id(go0, state_map, dfa_sets, q);
            int go1_id = get_state_id(go1, state_map, dfa_sets, q);

            tmp.push_back({cid, is_accept, go0_id, go1_id});
        }

        // 添加陷阱态
        if (state_map.find({}) == state_map.end()) {
            int tid = (int) tmp.size();
            state_map[{}] = tid;
            dfa_sets.push_back({});
            tmp.push_back({tid, false, tid, tid});
        }
        int trap_id = state_map[{}];
        for (auto &st: tmp) {
            if (st.t0 < 0) st.t0 = trap_id;
            if (st.t1 < 0) st.t1 = trap_id;
        }

        DFA dfa;
        dfa.states.resize(tmp.size());
        for (auto &st: tmp) {
            dfa.states[st.id].id = st.id;
            dfa.states[st.id].accept = st.acc;
            dfa.states[st.id].t0 = st.t0;
            dfa.states[st.id].t1 = st.t1;
        }
        dfa.start = 0;
        dfa.trap = trap_id;
        return dfa;
    }

private:
    const NFA &infa; ///< 输入NFA

    /**
     * @brief 获取子集对应的DFA状态ID，若无则创建新状态
     * @param S 当前NFA子集
     * @param m 子集到DFA状态ID映射
     * @param d dfa_sets集合
     * @param q 待处理队列
     * @return DFA状态ID
     */
    int get_state_id(const set<int> &S, map<set<int>, int> &m, vector<set<int> > &d, queue<set<int> > &q) {
        if (S.empty()) return -1;
        if (m.find(S) == m.end()) {
            int id = (int) d.size();
            m[S] = id;
            d.push_back(S);
            q.push(S);
            return id;
        }
        return m[S];
    }

    /**
     * @brief 根据输入字符c从NFA状态集合cur转移得到新集合
     * @param cur 当前NFA状态子集
     * @param c 输入字符('0'或'1')
     * @return 转移后的NFA状态子集
     */
    set<int> move_set(const set<int> &cur, char c) {
        set<int> res;
        for (auto s: cur) {
            auto it = infa.states[s].trans.find(c);
            if (it != infa.states[s].trans.end()) {
                for (auto nxt: it->second) {
                    res.insert(nxt);
                }
            }
        }
        return res;
    }
};

/**
 * @brief DFA最小化类。
 *
 * 使用Hopcroft或类似算法对DFA进行最小化，减少状态数。
 */
class DFAMinimizer {
public:
    /**
     * @brief 构造函数
     * @param d 待最小化的DFA
     */
    DFAMinimizer(const DFA &d) : idfa(d) {}

    /**
     * @brief 对DFA进行最小化
     * @return 最小化后的DFA
     */
    DFA minimize() {
        vector<int> partition(idfa.states.size());
        for (int i = 0; i < (int) idfa.states.size(); i++) {
            partition[i] = idfa.states[i].accept ? 1 : 0;
        }

        bool changed = true;
        while (changed) {
            changed = false;
            map<array<int, 3>, vector<int> > groups;
            for (int i = 0; i < (int) idfa.states.size(); i++) {
                int c0 = partition[idfa.states[i].t0];
                int c1 = partition[idfa.states[i].t1];
                array<int, 3> key = {partition[i], c0, c1};
                groups[key].push_back(i);
            }

            int new_class = 0;
            vector<int> new_partition(idfa.states.size());
            for (auto &g: groups) {
                for (auto s: g.second) {
                    new_partition[s] = new_class;
                }
                new_class++;
            }
            if (new_partition != partition) {
                partition = new_partition;
                changed = true;
            }
        }

        int count_classes = *max_element(partition.begin(), partition.end()) + 1;
        vector<int> repr(count_classes, -1);
        for (int i = 0; i < (int) partition.size(); i++) {
            if (repr[partition[i]] < 0) repr[partition[i]] = i;
        }

        DFA md;
        md.states.resize(count_classes);
        int start_class = partition[idfa.start];
        md.start = start_class;
        int trap_class = partition[idfa.trap];

        for (int c = 0; c < count_classes; c++) {
            int s = repr[c];
            bool class_accept = false;
            for (int i = 0; i < (int) idfa.states.size(); i++) {
                if (partition[i] == c && idfa.states[i].accept) {
                    class_accept = true;
                    break;
                }
            }

            md.states[c].id = c;
            md.states[c].accept = class_accept;
            md.states[c].t0 = partition[idfa.states[s].t0];
            md.states[c].t1 = partition[idfa.states[s].t1];
        }
        md.trap = trap_class;
        return md;
    }

private:
    const DFA &idfa; ///< 输入的DFA引用
};

/**
 * @brief 移除DFA中的陷阱态。
 *
 * 将陷阱态相关的转移标记为-1，陷阱态状态本身也设为无效。
 */
class TrapRemover {
public:
    /**
     * @brief 构造函数
     * @param d 待处理的DFA
     */
    TrapRemover(DFA &d) : dfa(d) {}

    /**
     * @brief 移除陷阱态，清除对陷阱态的转移
     */
    void remove_trap() {
        int t = dfa.trap;
        if (t == -1) return;
        for (auto &st: dfa.states) {
            if (st.t0 == t) st.t0 = -1;
            if (st.t1 == t) st.t1 = -1;
        }
        if (t >= 0 && t < (int) dfa.states.size()) {
            dfa.states[t].accept = false;
            dfa.states[t].t0 = -1;
            dfa.states[t].t1 = -1;
        }
        dfa.trap = -1;
    }

private:
    DFA &dfa; ///< DFA引用
};

/**
 * @brief 使用最终的最小化DFA进行分析的辅助类，主要用于演示和输出。
 */
class DFAPrinter {
public:
    /**
     * @brief 构造函数
     * @param d 最小化后的DFA引用
     */
    DFAPrinter(DFA &d) : idfa(d) {}

    /**
     * @brief 打印DFA和对应的右线性文法(RG)。
     */
    void print_and_convert_to_RG() {
        cout << "      0 1\n";
        vector<string> qname = name_states();
        auto order = get_order(qname);

        for (auto &pr: order) {
            int i = pr.second;
            bool start_mark = (i == idfa.start);
            bool accept_mark = idfa.states[i].accept;
            cout << (start_mark ? "(s)" : "") << (accept_mark ? "(e)" : "") << qname[i] << " ";
            print_transition(idfa.states[i].t0, qname);
            cout << " ";
            print_transition(idfa.states[i].t1, qname);
            cout << "\n";
        }

        cout << "\n";

        for (auto &p: order) {
            int stid = p.second;
            string lhs = qname[stid];
            vector<string> productions;

            process_symbol(stid, qname, productions);

            for (auto &prod: productions) {
                cout << prod << "\n";
            }
        }
    }

private:
    DFA &idfa; ///< DFA引用

    /**
     * @brief 为DFA状态命名，例如q0,q1,q2...
     * @return 每个状态对应的名称数组
     */
    vector<string> name_states() {
        vector<string> qname(idfa.states.size(), "");
        vector<bool> visited(idfa.states.size(), false);

        if (idfa.start >= 0 && idfa.start < (int) idfa.states.size()) {
            queue<int> Q;
            Q.push(idfa.start);
            visited[idfa.start] = true;
            qname[idfa.start] = "q0";
            int qid_count = 1;

            while (!Q.empty()) {
                int u = Q.front();
                Q.pop();
                int nxts[2] = {idfa.states[u].t0, idfa.states[u].t1};
                for (int i = 0; i < 2; i++) {
                    int v = nxts[i];
                    if (v >= 0 && v < (int) idfa.states.size() && !visited[v]) {
                        visited[v] = true;
                        qname[v] = "q" + to_string(qid_count++);
                        Q.push(v);
                    }
                }
            }
        }
        return qname;
    }

    /**
     * @brief 根据qX名称对状态排序
     * @param qname 状态名称数组
     * @return 排序后的(pair<索引,状态ID>)列表
     */
    vector<pair<int, int> > get_order(const vector<string> &qname) {
        vector<pair<int, int> > order;
        for (int i = 0; i < (int) qname.size(); i++) {
            if (!qname[i].empty() && qname[i][0] == 'q') {
                int idx = stoi(qname[i].substr(1));
                order.push_back({idx, i});
            }
        }
        sort(order.begin(), order.end());
        return order;
    }

    /**
     * @brief 打印状态转移，如果无效则打印N
     * @param tid 转移状态ID
     * @param qname 状态名称数组
     */
    void print_transition(int tid, const vector<string> &qname) {
        if (tid == -1 || tid >= (int) qname.size() || qname[tid].empty()) {
            cout << "N";
        } else {
            cout << qname[tid];
        }
    }

    /**
     * @brief 为给定状态生成右线性文法的产生式
     * @param stid 状态ID
     * @param qname 状态名称数组
     * @param productions 收集产生式的向量
     */
    void process_symbol(int stid, const vector<string> &qname, vector<string> &productions) {
        vector<string> production0;
        for (int c = 0; c < 2; c++) {
            int nxt = (c == 0) ? idfa.states[stid].t0 : idfa.states[stid].t1;
            if (nxt == -1 || nxt >= (int) qname.size() || qname[nxt].empty()) {
                continue;
            }

            string lhs = qname[stid];
            bool next_is_accept = idfa.states[nxt].accept;

            if (!(next_is_accept && idfa.states[nxt].t0 == -1 && idfa.states[nxt].t1 == -1)) {
                productions.push_back(lhs + "->" + to_string(c) + qname[nxt]);
            }
            if (next_is_accept) {
                production0.push_back(lhs + "->" + to_string(c));
            }
        }
        for (auto &p: production0) {
            productions.push_back(p);
        }
    }
};

/**
 * @brief 使用DFA来判断给定字符串是否被该正则表示的语言接受。
 * @param dfa 已构造并最小化的DFA
 * @param input 待判定的字符串，仅包含字符'0'和'1'
 * @return 若输入字符串被DFA接受则返回true，否则返回false
 */
bool check_membership(const DFA &dfa, const string &input) {
    int current_state = dfa.start;
    for (char c: input) {
        int next_state;
        if (c == '0') {
            next_state = dfa.states[current_state].t0;
        } else if (c == '1') {
            next_state = dfa.states[current_state].t1;
        } else {
            // 非法输入字符
            return false;
        }
        if (next_state == -1) {
            return false;
        }
        current_state = next_state;
    }
    return dfa.states[current_state].accept;
}

/**
 * @brief main函数入口。
 *
 * 流程：
 * 1. 从输入读取一个正则表达式字符串，如"(0+1)*0"
 * 2. 解析正则表达式得到AST
 * 3. 使用Thompson构造得ε-NFA，再进行ε消除得到NFA
 * 4. 使用子集构造法将NFA转为DFA
 * 5. 最小化DFA并移除陷阱态
 * 6. 从输入读取测试字符串并判定是否被DFA接受
 */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string re;
    cin >> re;

    RegexParser parser(re);
    RegexNode *root = parser.parse();

    Thompson th(root);
    NFA enfa = th.build();

    EpsilonRemover er(enfa);
    NFA nfa = er.remove();

    SubsetConstruction sc(nfa);
    DFA dfa = sc.convert();

    DFAMinimizer dm(dfa);
    DFA mdfa = dm.minimize();

    TrapRemover tr(mdfa);
    tr.remove_trap();

    DFAPrinter dfa_printer(mdfa);
    dfa_printer.print_and_convert_to_RG();

    string test_str;
    cin >> test_str;
    if (check_membership(mdfa, test_str)) {
        cout << "true\n";
    } else {
        cout << "false\n";
    }

    return 0;
}
