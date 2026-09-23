#define _GNU_SOURCE
#include "command_place.h"
#include "global_basic.h"
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <time.h>

#undef OBJ_BITS
#define GID_BITS 20
#define OBJ_BITS 20
#define OBJ_MASK ((uint64_t)((1ULL << OBJ_BITS) - 1ULL))
#define GID_MASK ((uint64_t)((1ULL << GID_BITS) - 1ULL))
#define STACK_WORDS 16

typedef struct __attribute__((packed)) {
    uint64_t ctxgid;
    uint32_t obj;
} sorted_rec_t;

typedef enum {
    RANK_WLS = 0,
    RANK_SPARSE = 1,
    RANK_SUM = 2
} rank_mode_t;

typedef enum {
    DISTANCE_FORMAT_AUTO = 0,
    DISTANCE_FORMAT_ANI = 1,
    DISTANCE_FORMAT_MATRIX = 2,
    DISTANCE_FORMAT_PHYLIP = 3
} distance_format_t;

typedef struct {
    char *name;
    double branch;
    int parent;
    int first_child;
    int next_sibling;
    int ref_id;
} node_t;

typedef struct {
    node_t *nodes;
    int n;
    int cap;
    int root;
} tree_t;

typedef struct {
    const char *s;
    size_t pos;
    size_t len;
    tree_t *tree;
} parser_t;

typedef struct {
    char **ids;
    char **samples;
    int n;
    int cap;
} refs_t;

typedef struct {
    char *name;
    double *dist;
    unsigned char *has_dist;
    int n_dist;
    int sketch_index;
} query_t;

typedef struct {
    query_t *items;
    int n;
    int cap;
    int preloaded;
} queries_t;

typedef struct {
    int to;
    int next;
    double len;
} adj_edge_t;

typedef struct {
    int edge_child;
    int edge_num;
    int clade_size;
    double rss;
    double pendant_length;
    double distal_from_child;
    double distal_from_parent;
    double edge_length;
    int wls_rank;
    double wls_delta;
    double wls_delta_ratio;
    double wls_rel_weight;
    double wls_top_weight;
    double wls_second_weight;
    double wls_entropy;
    double sparse_score;
    double sparse_score_per_ctx;
    int sparse_rank;
    int sparse_key_size;
    int sparse_uses_complement;
    int sparse_informative_ctx;
    int sparse_observed_ctx;
    uint64_t query_ctx;
    int sparse_used_ctx;
    double final_key;
} placement_t;

typedef struct {
    char *data;
    size_t len;
    size_t cap;
} strbuf_t;

typedef struct {
    tree_t *tree;
    refs_t *refs;
    queries_t *queries;
    int n_words;
    uint64_t *node_bits;
    double *node_dist;
    const char *weight_mode;
    rank_mode_t rank_mode;
    int top_n;
    int candidate_near;
    int sparse_nearest;
    int n_threads;
    double confidence_beta;
    double ambiguous_ratio;
    int min_sparse_informative;
    int sparse_enabled;
    sorted_rec_t *sorted;
    size_t n_sorted;
    uint64_t *qidx;
    size_t n_qidx;
    uint64_t *qcomb;
    size_t n_qcomb;
    char **blocks;
    placement_t **top_placements;
    int *n_top_placements;
    char **warnings;
} place_ctx_t;

typedef struct {
    place_ctx_t *ctx;
    int tid;
} worker_arg_t;

static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static char *xstrdup(const char *s) {
    char *out = strdup(s ? s : "");
    if (!out) die("strdup failed");
    return out;
}

static void *xcalloc(size_t n, size_t size) {
    void *p = calloc(n ? n : 1, size ? size : 1);
    if (!p) die("calloc failed for %zu x %zu", n, size);
    return p;
}

static void *xmalloc(size_t size) {
    void *p = malloc(size ? size : 1);
    if (!p) die("malloc failed for %zu bytes", size);
    return p;
}

static void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size ? size : 1);
    if (!p) die("realloc failed for %zu bytes", size);
    return p;
}

static char *read_text_file(const char *path, size_t *len_out) {
    struct stat st;
    if (stat(path, &st) != 0) die("stat failed for %s: %s", path, strerror(errno));
    FILE *fp = fopen(path, "rb");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    char *buf = (char *)xcalloc((size_t)st.st_size + 1, 1);
    size_t n = fread(buf, 1, (size_t)st.st_size, fp);
    if (n != (size_t)st.st_size) die("read failed for %s", path);
    fclose(fp);
    buf[n] = '\0';
    if (len_out) *len_out = n;
    return buf;
}

static void *read_binary_file(const char *path, size_t elem_size, size_t *n_out) {
    struct stat st;
    if (stat(path, &st) != 0) die("stat failed for %s: %s", path, strerror(errno));
    if ((size_t)st.st_size % elem_size != 0) {
        die("bad file size for %s: %jd is not divisible by %zu", path, (intmax_t)st.st_size, elem_size);
    }
    FILE *fp = fopen(path, "rb");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    size_t n = (size_t)st.st_size / elem_size;
    void *buf = xmalloc((size_t)st.st_size);
    if (n && fread(buf, elem_size, n, fp) != n) die("read failed for %s", path);
    fclose(fp);
    if (n_out) *n_out = n;
    return buf;
}

static char *path_join(const char *dir, const char *file) {
    size_t n = strlen(dir) + strlen(file) + 2;
    char *out = (char *)xmalloc(n);
    snprintf(out, n, "%s/%s", dir, file);
    return out;
}

typedef struct {
    dim_sketch_stat_t stat;
    char **names;
    uint64_t *index;
    size_t n_index;
    uint64_t *comb;
    size_t n_comb;
} query_sketch_part_t;

static void query_sketch_part_destroy(query_sketch_part_t *part) {
    if (!part) return;
    for (int i = 0; i < part->stat.infile_num; i++) free(part->names ? part->names[i] : NULL);
    free(part->names);
    free(part->index);
    free(part->comb);
    memset(part, 0, sizeof(*part));
}

static void check_query_sketch_stats(const dim_sketch_stat_t *expected,
                                     const dim_sketch_stat_t *actual,
                                     const char *dir) {
    if (expected->hash_id != actual->hash_id ||
        expected->koc != actual->koc ||
        expected->conflict != actual->conflict ||
        expected->coden_len != actual->coden_len ||
        expected->klen != actual->klen ||
        expected->hclen != actual->hclen ||
        expected->holen != actual->holen ||
        expected->drfold != actual->drfold) {
        die("query sketch %s is incompatible with the first --query-sketch-list entry", dir);
    }
}

static query_sketch_part_t load_query_sketch_part(const char *dir) {
    query_sketch_part_t part = {0};
    size_t stat_len = 0;
    char *stat_path = path_join(dir, "lcofiles.stat");
    char *stat_data = read_text_file(stat_path, &stat_len);
    if (stat_len < sizeof(part.stat)) {
        die("invalid lcofiles.stat in %s", dir);
    }
    memcpy(&part.stat, stat_data, sizeof(part.stat));
    if (part.stat.infile_num <= 0 ||
        (size_t)part.stat.infile_num > (SIZE_MAX - sizeof(part.stat)) / PATHLEN ||
        stat_len != sizeof(part.stat) + (size_t)part.stat.infile_num * PATHLEN) {
        die("invalid lcofiles.stat record count in %s", dir);
    }
    part.names = (char **)xcalloc((size_t)part.stat.infile_num, sizeof(*part.names));
    for (int i = 0; i < part.stat.infile_num; i++) {
        const char *stored_name = stat_data + sizeof(part.stat) + (size_t)i * PATHLEN;
        if (!memchr(stored_name, '\0', PATHLEN)) {
            die("unterminated sample name in lcofiles.stat: %s", dir);
        }
        part.names[i] = xstrdup(stored_name);
    }
    free(stat_data);
    free(stat_path);

    char *index_path = path_join(dir, "comblco.index");
    char *comb_path = path_join(dir, "comblco");
    part.index = (uint64_t *)read_binary_file(index_path, sizeof(uint64_t), &part.n_index);
    part.comb = (uint64_t *)read_binary_file(comb_path, sizeof(uint64_t), &part.n_comb);
    free(index_path);
    free(comb_path);
    if (part.n_index != (size_t)part.stat.infile_num + 1 || part.index[0] != 0 ||
        part.index[part.n_index - 1] != part.n_comb) {
        die("invalid comblco index in %s", dir);
    }
    for (size_t i = 1; i < part.n_index; i++) {
        if (part.index[i] < part.index[i - 1]) {
            die("non-monotonic comblco index in %s", dir);
        }
    }
    return part;
}

static char *trim_query_sketch_list_line(char *line) {
    while (*line && isspace((unsigned char)*line)) line++;
    char *end = line + strlen(line);
    while (end > line && isspace((unsigned char)end[-1])) *--end = '\0';
    return line;
}

static void load_query_sketch_list(const char *list_path,
                                   uint64_t **index_out, size_t *n_index_out,
                                   uint64_t **comb_out, size_t *n_comb_out,
                                   char ***names_out, size_t *n_names_out) {
    FILE *fp = fopen(list_path, "r");
    if (!fp) die("open failed for %s: %s", list_path, strerror(errno));
    query_sketch_part_t *parts = NULL;
    size_t n_parts = 0, cap_parts = 0, total_samples = 0, total_comb = 0;
    dim_sketch_stat_t expected = {0};
    char *line = NULL;
    size_t line_cap = 0;
    ssize_t line_len;
    while ((line_len = getline(&line, &line_cap, fp)) >= 0) {
        (void)line_len;
        char *dir = trim_query_sketch_list_line(line);
        if (!*dir || *dir == '#') continue;
        if (n_parts == cap_parts) {
            cap_parts = cap_parts ? cap_parts * 2 : 16;
            parts = (query_sketch_part_t *)xrealloc(parts, cap_parts * sizeof(*parts));
        }
        parts[n_parts] = load_query_sketch_part(dir);
        if (n_parts == 0) expected = parts[n_parts].stat;
        else check_query_sketch_stats(&expected, &parts[n_parts].stat, dir);
        if ((size_t)parts[n_parts].stat.infile_num > SIZE_MAX - total_samples ||
            parts[n_parts].n_comb > SIZE_MAX - total_comb) {
            die("query sketch list is too large");
        }
        total_samples += (size_t)parts[n_parts].stat.infile_num;
        total_comb += parts[n_parts].n_comb;
        n_parts++;
    }
    free(line);
    fclose(fp);
    if (n_parts == 0) die("no query sketch directories listed in %s", list_path);
    if (total_samples == SIZE_MAX || total_comb > SIZE_MAX / sizeof(uint64_t)) {
        die("query sketch list is too large");
    }

    uint64_t *all_index = (uint64_t *)xcalloc(total_samples + 1, sizeof(*all_index));
    uint64_t *all_comb = (uint64_t *)xmalloc(total_comb * sizeof(*all_comb));
    char **all_names = (char **)xcalloc(total_samples, sizeof(*all_names));
    size_t sample_offset = 0, comb_offset = 0;
    for (size_t p = 0; p < n_parts; p++) {
        query_sketch_part_t *part = &parts[p];
        size_t n_samples = (size_t)part->stat.infile_num;
        for (size_t i = 1; i <= n_samples; i++) {
            all_index[sample_offset + i] = comb_offset + part->index[i];
        }
        for (size_t i = 0; i < n_samples; i++) {
            all_names[sample_offset + i] = part->names[i];
            part->names[i] = NULL;
        }
        if (part->n_comb) {
            memcpy(all_comb + comb_offset, part->comb, part->n_comb * sizeof(*all_comb));
        }
        sample_offset += n_samples;
        comb_offset += part->n_comb;
        query_sketch_part_destroy(part);
    }
    free(parts);
    *index_out = all_index;
    *n_index_out = total_samples + 1;
    *comb_out = all_comb;
    *n_comb_out = total_comb;
    *names_out = all_names;
    *n_names_out = total_samples;
}

static void load_query_sketch_dir(const char *dir,
                                  uint64_t **index_out, size_t *n_index_out,
                                  uint64_t **comb_out, size_t *n_comb_out,
                                  char ***names_out, size_t *n_names_out) {
    query_sketch_part_t part = load_query_sketch_part(dir);
    *index_out = part.index;
    *n_index_out = part.n_index;
    *comb_out = part.comb;
    *n_comb_out = part.n_comb;
    *names_out = part.names;
    *n_names_out = (size_t)part.stat.infile_num;
    part.index = NULL;
    part.comb = NULL;
    part.names = NULL;
    query_sketch_part_destroy(&part);
}

static void free_query_sketch_names(char **names, size_t n_names) {
    if (!names) return;
    for (size_t i = 0; i < n_names; i++) free(names[i]);
    free(names);
}

static void strbuf_appendf(strbuf_t *b, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    va_list ap2;
    va_copy(ap2, ap);
    int need = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (need < 0) die("vsnprintf failed");
    size_t add = (size_t)need;
    if (b->len + add + 1 > b->cap) {
        size_t cap = b->cap ? b->cap : 4096;
        while (b->len + add + 1 > cap) cap *= 2;
        b->data = (char *)xrealloc(b->data, cap);
        b->cap = cap;
    }
    vsnprintf(b->data + b->len, b->cap - b->len, fmt, ap2);
    va_end(ap2);
    b->len += add;
}

static void tree_add_capacity(tree_t *t) {
    if (t->n == t->cap) {
        t->cap = t->cap ? t->cap * 2 : 512;
        t->nodes = (node_t *)xrealloc(t->nodes, (size_t)t->cap * sizeof(node_t));
    }
}

static int tree_new_node(tree_t *t, int parent) {
    tree_add_capacity(t);
    int id = t->n++;
    memset(&t->nodes[id], 0, sizeof(t->nodes[id]));
    t->nodes[id].parent = parent;
    t->nodes[id].first_child = -1;
    t->nodes[id].next_sibling = -1;
    t->nodes[id].ref_id = -1;
    return id;
}

static void tree_append_child(tree_t *t, int parent, int child) {
    if (t->nodes[parent].first_child < 0) {
        t->nodes[parent].first_child = child;
        return;
    }
    int x = t->nodes[parent].first_child;
    while (t->nodes[x].next_sibling >= 0) x = t->nodes[x].next_sibling;
    t->nodes[x].next_sibling = child;
}

static void skip_ws(parser_t *p) {
    while (p->pos < p->len && isspace((unsigned char)p->s[p->pos])) p->pos++;
}

static char *parse_label(parser_t *p) {
    skip_ws(p);
    if (p->pos >= p->len) return xstrdup("");
    if (p->s[p->pos] == '\'') {
        p->pos++;
        size_t start = p->pos;
        while (p->pos < p->len && p->s[p->pos] != '\'') p->pos++;
        size_t n = p->pos - start;
        char *out = (char *)xcalloc(n + 1, 1);
        memcpy(out, p->s + start, n);
        if (p->pos < p->len && p->s[p->pos] == '\'') p->pos++;
        return out;
    }
    size_t start = p->pos;
    while (p->pos < p->len) {
        char c = p->s[p->pos];
        if (c == ':' || c == ',' || c == ')' || c == '(' || c == ';' || isspace((unsigned char)c)) break;
        p->pos++;
    }
    size_t n = p->pos - start;
    char *out = (char *)xcalloc(n + 1, 1);
    memcpy(out, p->s + start, n);
    return out;
}

static int parse_subtree(parser_t *p, int parent) {
    skip_ws(p);
    int node = tree_new_node(p->tree, parent);
    if (p->pos < p->len && p->s[p->pos] == '(') {
        p->pos++;
        for (;;) {
            int child = parse_subtree(p, node);
            tree_append_child(p->tree, node, child);
            skip_ws(p);
            if (p->pos >= p->len) die("unexpected end while parsing Newick");
            if (p->s[p->pos] == ',') {
                p->pos++;
                continue;
            }
            if (p->s[p->pos] == ')') {
                p->pos++;
                break;
            }
            die("unexpected Newick character '%c'", p->s[p->pos]);
        }
        char *label = parse_label(p);
        if (label[0]) p->tree->nodes[node].name = label;
        else free(label);
    } else {
        char *label = parse_label(p);
        p->tree->nodes[node].name = label;
    }
    skip_ws(p);
    if (p->pos < p->len && p->s[p->pos] == ':') {
        p->pos++;
        char *end = NULL;
        errno = 0;
        double v = strtod(p->s + p->pos, &end);
        if (errno || end == p->s + p->pos) die("bad branch length near byte %zu", p->pos);
        p->tree->nodes[node].branch = v;
        p->pos = (size_t)(end - p->s);
    }
    return node;
}

static tree_t read_newick(const char *path) {
    size_t len = 0;
    char *text = read_text_file(path, &len);
    tree_t tree = {0};
    parser_t parser = {.s = text, .pos = 0, .len = len, .tree = &tree};
    tree.root = parse_subtree(&parser, -1);
    skip_ws(&parser);
    if (parser.pos < parser.len && parser.s[parser.pos] == ';') parser.pos++;
    free(text);
    return tree;
}

static void refs_push(refs_t *refs, const char *id, const char *sample) {
    if (refs->n == refs->cap) {
        refs->cap = refs->cap ? refs->cap * 2 : 256;
        refs->ids = (char **)xrealloc(refs->ids, (size_t)refs->cap * sizeof(char *));
        refs->samples = (char **)xrealloc(refs->samples, (size_t)refs->cap * sizeof(char *));
    }
    refs->ids[refs->n] = xstrdup(id);
    refs->samples[refs->n] = xstrdup(sample);
    refs->n++;
}

static refs_t read_idmap(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    refs_t refs = {0};
    char *line = NULL;
    size_t len = 0;
    if (getline(&line, &len, fp) < 0) die("empty idmap: %s", path);
    while (getline(&line, &len, fp) >= 0) {
        line[strcspn(line, "\r\n")] = '\0';
        if (!line[0]) continue;
        char *tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = '\0';
        refs_push(&refs, line, tab + 1);
    }
    free(line);
    fclose(fp);
    return refs;
}

static int find_ref_by_id(const refs_t *refs, const char *id) {
    for (int i = 0; i < refs->n; i++) {
        if (strcmp(refs->ids[i], id) == 0) return i;
    }
    return -1;
}

static int find_ref_by_sample(const refs_t *refs, const char *sample) {
    for (int i = 0; i < refs->n; i++) {
        if (strcmp(refs->samples[i], sample) == 0) return i;
    }
    return -1;
}

static int find_ref_by_name(const refs_t *refs, const char *name) {
    int id = find_ref_by_sample(refs, name);
    if (id >= 0) return id;
    return find_ref_by_id(refs, name);
}

static void map_tree_refs(tree_t *tree, const refs_t *refs) {
    for (int i = 0; i < tree->n; i++) {
        if (tree->nodes[i].name) tree->nodes[i].ref_id = find_ref_by_id(refs, tree->nodes[i].name);
    }
}

static query_t *queries_find(queries_t *queries, const char *name) {
    for (int i = 0; i < queries->n; i++) {
        if (strcmp(queries->items[i].name, name) == 0) return &queries->items[i];
    }
    return NULL;
}

static query_t *queries_add(queries_t *queries, const char *name, int n_refs) {
    if (queries->n == queries->cap) {
        queries->cap = queries->cap ? queries->cap * 2 : 64;
        queries->items = (query_t *)xrealloc(queries->items, (size_t)queries->cap * sizeof(query_t));
    }
    query_t *q = &queries->items[queries->n++];
    memset(q, 0, sizeof(*q));
    q->name = xstrdup(name);
    q->dist = (double *)xcalloc((size_t)n_refs, sizeof(double));
    q->has_dist = (unsigned char *)xcalloc((size_t)n_refs, 1);
    q->sketch_index = -1;
    return q;
}

static query_t *queries_get_or_add(queries_t *queries, const char *name, int n_refs) {
    query_t *q = queries_find(queries, name);
    if (q) return q;
    if (queries->preloaded) return NULL;
    return queries_add(queries, name, n_refs);
}

static queries_t read_query_paths(const char *path, int n_refs) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    queries_t queries = {.preloaded = 1};
    char *line = NULL;
    size_t len = 0;
    while (getline(&line, &len, fp) >= 0) {
        line[strcspn(line, "\r\n")] = '\0';
        if (!line[0]) continue;
        queries_add(&queries, line, n_refs);
    }
    free(line);
    fclose(fp);
    return queries;
}

static queries_t read_query_names(char *const *names, size_t n_names, int n_refs) {
    queries_t queries = {.preloaded = 1};
    for (size_t i = 0; i < n_names; i++) {
        if (!names[i] || !names[i][0]) die("empty sample name in query sketch metadata");
        if (queries_find(&queries, names[i])) {
            die("duplicate sample name in query sketch metadata: %s", names[i]);
        }
        query_t *q = queries_add(&queries, names[i], n_refs);
        q->sketch_index = (int)i;
    }
    return queries;
}

static void validate_query_paths(const char *path, char *const *names, size_t n_names) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    char *line = NULL;
    size_t len = 0;
    size_t i = 0;
    while (getline(&line, &len, fp) >= 0) {
        line[strcspn(line, "\r\n")] = '\0';
        if (!line[0]) continue;
        if (i >= n_names) {
            die("--query-paths has more labels than query sketch metadata");
        }
        if (strcmp(line, names[i]) != 0) {
            die("--query-paths label %zu (%s) does not match query sketch sample (%s); "
                "query sketch metadata is authoritative",
                i + 1, line, names[i]);
        }
        i++;
    }
    free(line);
    fclose(fp);
    if (i != n_names) {
        die("--query-paths has %zu labels; query sketch metadata has %zu samples", i, n_names);
    }
}

static queries_t init_queries(const char *query_paths, char *const *query_names,
                              size_t n_query_names, int n_refs) {
    if (query_names) {
        if (query_paths) validate_query_paths(query_paths, query_names, n_query_names);
        return read_query_names(query_names, n_query_names, n_refs);
    }
    if (query_paths) return read_query_paths(query_paths, n_refs);
    return (queries_t){0};
}

static int split_tab(char *line, char **cols, int max_cols) {
    int n = 0;
    char *save = NULL;
    char *tok = strtok_r(line, "\t\r\n", &save);
    while (tok && n < max_cols) {
        cols[n++] = tok;
        tok = strtok_r(NULL, "\t\r\n", &save);
    }
    return n;
}

static int split_tab_keep_empty(char *line, char **cols, int max_cols) {
    line[strcspn(line, "\r\n")] = '\0';
    int n = 0;
    char *p = line;
    if (n < max_cols) cols[n++] = p;
    while (*p) {
        if (*p == '\t') {
            *p = '\0';
            if (n < max_cols) cols[n++] = p + 1;
        }
        p++;
    }
    return n;
}

static int split_ws(char *line, char **cols, int max_cols) {
    int n = 0;
    char *save = NULL;
    char *tok = strtok_r(line, " \t\r\n", &save);
    while (tok && n < max_cols) {
        cols[n++] = tok;
        tok = strtok_r(NULL, " \t\r\n", &save);
    }
    return n;
}

static int header_index(char **cols, int n, const char *name) {
    for (int i = 0; i < n; i++) {
        if (strcmp(cols[i], name) == 0) return i;
    }
    return -1;
}

static distance_format_t parse_distance_format(const char *s) {
    if (!s || strcmp(s, "auto") == 0) return DISTANCE_FORMAT_AUTO;
    if (strcmp(s, "ani") == 0 || strcmp(s, "detail") == 0 || strcmp(s, "long") == 0)
        return DISTANCE_FORMAT_ANI;
    if (strcmp(s, "matrix") == 0 || strcmp(s, "tsv-matrix") == 0 || strcmp(s, "tsv") == 0)
        return DISTANCE_FORMAT_MATRIX;
    if (strcmp(s, "phylip") == 0)
        return DISTANCE_FORMAT_PHYLIP;
    die("unknown --distance-format '%s' (expected auto, ani, matrix, or phylip)", s);
    return DISTANCE_FORMAT_AUTO;
}

static int is_integer_line(const char *line) {
    const unsigned char *p = (const unsigned char *)line;
    while (*p && isspace(*p)) p++;
    if (!isdigit(*p)) return 0;
    while (*p && isdigit(*p)) p++;
    while (*p && isspace(*p)) p++;
    return *p == '\0';
}

static distance_format_t detect_distance_format(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    char *line = NULL;
    size_t len = 0;
    if (getline(&line, &len, fp) < 0) die("empty distance table: %s", path);
    fclose(fp);

    char *header = xstrdup(line);
    char *hcols[128] = {0};
    int hn = split_tab(header, hcols, 128);
    const int has_long_cols = header_index(hcols, hn, "Qry") >= 0 &&
                              header_index(hcols, hn, "Ref") >= 0 &&
                              header_index(hcols, hn, "Distance") >= 0;
    free(header);

    distance_format_t fmt = DISTANCE_FORMAT_MATRIX;
    if (has_long_cols) {
        fmt = DISTANCE_FORMAT_ANI;
    } else if (is_integer_line(line)) {
        fmt = DISTANCE_FORMAT_PHYLIP;
    } else if (strchr(line, '\t')) {
        fmt = DISTANCE_FORMAT_MATRIX;
    } else {
        die("cannot auto-detect distance format for %s; use --distance-format ani|matrix|phylip", path);
    }
    free(line);
    return fmt;
}

static queries_t read_distances_ani(const char *path, const refs_t *refs, const char *query_paths,
                                    char *const *query_names, size_t n_query_names) {
    queries_t queries = init_queries(query_paths, query_names, n_query_names, refs->n);
    FILE *fp = fopen(path, "r");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    char *line = NULL;
    size_t len = 0;
    if (getline(&line, &len, fp) < 0) die("empty distance table: %s", path);
    char *header = xstrdup(line);
    char *hcols[128] = {0};
    int hn = split_tab(header, hcols, 128);
    int qi = header_index(hcols, hn, "Qry");
    int ri = header_index(hcols, hn, "Ref");
    int di = header_index(hcols, hn, "Distance");
    if (qi < 0 || ri < 0 || di < 0) {
        die("distance table must contain Qry, Ref, and Distance columns: %s", path);
    }
    while (getline(&line, &len, fp) >= 0) {
        char *cols[128] = {0};
        int n = split_tab(line, cols, 128);
        if (n <= qi || n <= ri || n <= di) continue;
        int ref_id = find_ref_by_name(refs, cols[ri]);
        if (ref_id < 0) continue;
        query_t *q = queries_get_or_add(&queries, cols[qi], refs->n);
        if (!q) continue;
        if (!q->has_dist[ref_id]) q->n_dist++;
        q->dist[ref_id] = atof(cols[di]);
        q->has_dist[ref_id] = 1;
    }
    free(header);
    free(line);
    fclose(fp);
    return queries;
}

static queries_t read_distances_matrix(const char *path, const refs_t *refs, const char *query_paths,
                                       char *const *query_names, size_t n_query_names) {
    queries_t queries = init_queries(query_paths, query_names, n_query_names, refs->n);
    FILE *fp = fopen(path, "r");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    char *line = NULL;
    size_t len = 0;
    if (getline(&line, &len, fp) < 0) die("empty distance matrix: %s", path);

    char *header = xstrdup(line);
    char **hcols = (char **)xcalloc((size_t)refs->n + 2, sizeof(char *));
    int hn = split_tab_keep_empty(header, hcols, refs->n + 2);
    if (hn < 2) die("distance matrix header must contain at least one reference column: %s", path);

    int *col_ref = (int *)xmalloc((size_t)hn * sizeof(int));
    for (int i = 0; i < hn; i++) col_ref[i] = -1;
    int mapped_refs = 0;
    for (int i = 1; i < hn; i++) {
        if (!hcols[i][0]) continue;
        col_ref[i] = find_ref_by_name(refs, hcols[i]);
        if (col_ref[i] >= 0) mapped_refs++;
    }
    if (mapped_refs == 0)
        die("distance matrix header does not contain any idmap reference names or ids: %s", path);

    char **cols = (char **)xcalloc((size_t)hn + 1, sizeof(char *));
    while (getline(&line, &len, fp) >= 0) {
        int n = split_tab_keep_empty(line, cols, hn + 1);
        if (n <= 1 || !cols[0][0]) continue;
        query_t *q = queries_get_or_add(&queries, cols[0], refs->n);
        if (!q) continue;
        for (int i = 1; i < hn && i < n; i++) {
            if (col_ref[i] < 0 || !cols[i][0]) continue;
            char *end = NULL;
            errno = 0;
            double d = strtod(cols[i], &end);
            if (errno || end == cols[i] || !isfinite(d)) continue;
            if (!q->has_dist[col_ref[i]]) q->n_dist++;
            q->dist[col_ref[i]] = d;
            q->has_dist[col_ref[i]] = 1;
        }
    }

    free(cols);
    free(col_ref);
    free(hcols);
    free(header);
    free(line);
    fclose(fp);
    return queries;
}

static queries_t read_distances_phylip(const char *path, const refs_t *refs, const char *query_paths,
                                       char *const *query_names, size_t n_query_names) {
    queries_t queries = init_queries(query_paths, query_names, n_query_names, refs->n);
    FILE *fp = fopen(path, "r");
    if (!fp) die("open failed for %s: %s", path, strerror(errno));
    char *line = NULL;
    size_t len = 0;
    if (getline(&line, &len, fp) < 0) die("empty PHYLIP distance matrix: %s", path);
    char *end = NULL;
    errno = 0;
    long n_long = strtol(line, &end, 10);
    if (errno || n_long <= 0 || n_long > 1000000)
        die("bad PHYLIP matrix size in %s", path);
    int n = (int)n_long;

    char **labels = (char **)xcalloc((size_t)n, sizeof(char *));
    double *values = (double *)xmalloc((size_t)n * (size_t)n * sizeof(double));
    char **cols = (char **)xcalloc((size_t)n + 2, sizeof(char *));
    for (int row = 0; row < n; row++) {
        if (getline(&line, &len, fp) < 0)
            die("PHYLIP matrix ended early at row %d/%d: %s", row + 1, n, path);
        int nc = split_ws(line, cols, n + 2);
        if (nc < n + 1)
            die("PHYLIP row %d has %d columns; expected label plus %d distances: %s",
                row + 1, nc, n, path);
        labels[row] = xstrdup(cols[0]);
        for (int col = 0; col < n; col++) {
            char *v_end = NULL;
            errno = 0;
            double d = strtod(cols[col + 1], &v_end);
            if (errno || v_end == cols[col + 1] || !isfinite(d))
                die("bad PHYLIP distance at row %d col %d in %s", row + 1, col + 1, path);
            values[(size_t)row * (size_t)n + (size_t)col] = d;
        }
    }

    int *col_ref = (int *)xmalloc((size_t)n * sizeof(int));
    int mapped_refs = 0;
    for (int col = 0; col < n; col++) {
        col_ref[col] = find_ref_by_name(refs, labels[col]);
        if (col_ref[col] >= 0) mapped_refs++;
    }
    if (mapped_refs == 0)
        die("PHYLIP matrix labels do not contain any idmap reference names or ids: %s", path);

    for (int row = 0; row < n; row++) {
        if (find_ref_by_name(refs, labels[row]) >= 0 && !query_paths)
            continue;
        query_t *q = queries_get_or_add(&queries, labels[row], refs->n);
        if (!q) continue;
        for (int col = 0; col < n; col++) {
            if (col_ref[col] < 0) continue;
            if (!q->has_dist[col_ref[col]]) q->n_dist++;
            q->dist[col_ref[col]] = values[(size_t)row * (size_t)n + (size_t)col];
            q->has_dist[col_ref[col]] = 1;
        }
    }

    for (int i = 0; i < n; i++) free(labels[i]);
    free(labels);
    free(values);
    free(cols);
    free(col_ref);
    free(line);
    fclose(fp);
    return queries;
}

static queries_t read_distances(const char *path, const refs_t *refs, const char *query_paths,
                                char *const *query_names, size_t n_query_names,
                                distance_format_t format) {
    if (format == DISTANCE_FORMAT_AUTO)
        format = detect_distance_format(path);
    switch (format) {
    case DISTANCE_FORMAT_ANI:
        return read_distances_ani(path, refs, query_paths, query_names, n_query_names);
    case DISTANCE_FORMAT_MATRIX:
        return read_distances_matrix(path, refs, query_paths, query_names, n_query_names);
    case DISTANCE_FORMAT_PHYLIP:
        return read_distances_phylip(path, refs, query_paths, query_names, n_query_names);
    case DISTANCE_FORMAT_AUTO:
    default:
        break;
    }
    die("unsupported distance format for %s", path);
    return (queries_t){0};
}

static int popcount_words(const uint64_t *words, int n_words) {
    int total = 0;
    for (int i = 0; i < n_words; i++) total += __builtin_popcountll(words[i]);
    return total;
}

static int popcount_and(const uint64_t *a, const uint64_t *b, int n_words) {
    int total = 0;
    for (int i = 0; i < n_words; i++) total += __builtin_popcountll(a[i] & b[i]);
    return total;
}

static int bit_is_set(const uint64_t *words, int bit) {
    return (int)((words[bit / 64] >> (bit % 64)) & 1ULL);
}

static void bit_set(uint64_t *words, int bit) {
    words[bit / 64] |= 1ULL << (bit % 64);
}

static void fill_node_bits_rec(const tree_t *tree, int node, uint64_t *bits, int n_words) {
    uint64_t *my = bits + (size_t)node * n_words;
    int ref_id = tree->nodes[node].ref_id;
    if (ref_id >= 0) bit_set(my, ref_id);
    for (int child = tree->nodes[node].first_child; child >= 0; child = tree->nodes[child].next_sibling) {
        fill_node_bits_rec(tree, child, bits, n_words);
        uint64_t *cb = bits + (size_t)child * n_words;
        for (int w = 0; w < n_words; w++) my[w] |= cb[w];
    }
}

static uint64_t *build_node_bits(const tree_t *tree, int n_refs, int *n_words_out) {
    int n_words = (n_refs + 63) / 64;
    uint64_t *bits = (uint64_t *)xcalloc((size_t)tree->n * n_words, sizeof(uint64_t));
    fill_node_bits_rec(tree, tree->root, bits, n_words);
    *n_words_out = n_words;
    return bits;
}

static double nonnegative_branch(double x) {
    return x > 0.0 ? x : 0.0;
}

static int count_negative_branches(const tree_t *tree) {
    int n = 0;
    for (int i = 0; i < tree->n; i++) {
        if (i != tree->root && tree->nodes[i].branch < 0.0) n++;
    }
    return n;
}

static void add_adj(adj_edge_t *edges, int *head, int *ep, int a, int b, double len) {
    edges[*ep].to = b;
    edges[*ep].len = len;
    edges[*ep].next = head[a];
    head[a] = (*ep)++;
}

static void build_adj(const tree_t *tree, int **head_out, adj_edge_t **edges_out) {
    int *head = (int *)xcalloc((size_t)tree->n, sizeof(int));
    for (int i = 0; i < tree->n; i++) head[i] = -1;
    adj_edge_t *edges = (adj_edge_t *)xcalloc((size_t)(2 * (tree->n - 1) + 2), sizeof(adj_edge_t));
    int ep = 0;
    for (int child = 0; child < tree->n; child++) {
        int parent = tree->nodes[child].parent;
        if (parent < 0) continue;
        double len = nonnegative_branch(tree->nodes[child].branch);
        add_adj(edges, head, &ep, parent, child, len);
        add_adj(edges, head, &ep, child, parent, len);
    }
    *head_out = head;
    *edges_out = edges;
}

static double *precompute_node_ref_dist(const tree_t *tree, const refs_t *refs) {
    int *head = NULL;
    adj_edge_t *edges = NULL;
    build_adj(tree, &head, &edges);
    double *dist = (double *)xcalloc((size_t)tree->n * refs->n, sizeof(double));
    int *stack = (int *)xcalloc((size_t)tree->n, sizeof(int));
    int *parent = (int *)xcalloc((size_t)tree->n, sizeof(int));
    double *sdist = (double *)xcalloc((size_t)tree->n, sizeof(double));
    for (int start = 0; start < tree->n; start++) {
        int sp = 0;
        stack[sp] = start;
        parent[sp] = -1;
        sdist[sp] = 0.0;
        sp++;
        while (sp) {
            sp--;
            int node = stack[sp];
            int par = parent[sp];
            double d = sdist[sp];
            int ref_id = tree->nodes[node].ref_id;
            if (ref_id >= 0) dist[(size_t)start * refs->n + ref_id] = d;
            for (int e = head[node]; e >= 0; e = edges[e].next) {
                int to = edges[e].to;
                if (to == par) continue;
                stack[sp] = to;
                parent[sp] = node;
                sdist[sp] = d + edges[e].len;
                sp++;
            }
        }
    }
    free(head);
    free(edges);
    free(stack);
    free(parent);
    free(sdist);
    return dist;
}

static double weight_for(double d, const char *mode) {
    const double eps = 1e-6;
    double x = d > eps ? d : eps;
    if (strcmp(mode, "equal") == 0) return 1.0;
    if (strcmp(mode, "inv_d") == 0) return 1.0 / x;
    if (strcmp(mode, "inv_d2") == 0) return 1.0 / (x * x);
    if (strcmp(mode, "inv_d_sqrt") == 0) return 1.0 / sqrt(x);
    if (strcmp(mode, "inv_d_cap") == 0) {
        double w = 1.0 / (x * x);
        return w < 1e6 ? w : 1e6;
    }
    if (strcmp(mode, "close_exp") == 0) return exp(-d / 0.005);
    die("unknown --weight-mode: %s", mode);
    return 1.0;
}

static void append_warning(char *warn, size_t n, const char *token) {
    if (!token || !token[0] || n == 0) return;
    size_t used = strlen(warn);
    if (used > 0) {
        if (used + 1 >= n) return;
        warn[used++] = ';';
        warn[used] = '\0';
    }
    if (used < n) snprintf(warn + used, n - used, "%s", token);
}

static void compute_wls_confidence(const place_ctx_t *ctx, placement_t *placements, int n_place) {
    if (n_place <= 0) return;
    double best = placements[0].rss;
    for (int i = 1; i < n_place; i++) {
        if (placements[i].rss < best) best = placements[i].rss;
    }
    double scale = fabs(best);
    if (scale < 1e-12) scale = 1e-12;

    double *weights = (double *)xcalloc((size_t)n_place, sizeof(double));
    double sum = 0.0;
    for (int i = 0; i < n_place; i++) {
        double delta = placements[i].rss - best;
        if (delta < 0.0) delta = 0.0;
        double ratio = delta / scale;
        double logw = -ctx->confidence_beta * ratio;
        double w = logw < -745.0 ? 0.0 : exp(logw);
        placements[i].wls_delta = delta;
        placements[i].wls_delta_ratio = ratio;
        weights[i] = w;
        sum += w;
    }

    double top1 = 0.0;
    double top2 = 0.0;
    double entropy = 0.0;
    for (int i = 0; i < n_place; i++) {
        double p = sum > 0.0 ? weights[i] / sum : 0.0;
        placements[i].wls_rel_weight = p;
        if (p > top1) {
            top2 = top1;
            top1 = p;
        } else if (p > top2) {
            top2 = p;
        }
        if (p > 0.0) entropy -= p * log(p);
    }
    if (n_place > 1) entropy /= log((double)n_place);
    for (int i = 0; i < n_place; i++) {
        placements[i].wls_top_weight = top1;
        placements[i].wls_second_weight = top2;
        placements[i].wls_entropy = entropy;
    }
    free(weights);
}

static void fit_edge(const double *y, const double *s, const double *w, int n, double edge_len,
                     double *rss_out, double *pendant_out, double *distal_child_out) {
    double sw = 0.0, sy = 0.0, ss = 0.0, ssy = 0.0, sss = 0.0;
    for (int i = 0; i < n; i++) {
        sw += w[i];
        sy += w[i] * y[i];
        ss += w[i] * s[i];
        ssy += w[i] * s[i] * y[i];
        sss += w[i] * s[i] * s[i];
    }
    double best_rss = DBL_MAX;
    double best_b = 0.0;
    double best_t = 0.0;
    double cand_b[8];
    double cand_t[8];
    int nc = 0;
    double det = sw * sss - ss * ss;
    if (det > 0.0) {
        double b = (sy * sss - ss * ssy) / det;
        double t = (sw * ssy - ss * sy) / det;
        if (b < 0.0) b = 0.0;
        if (t < 0.0) t = 0.0;
        if (t > edge_len) t = edge_len;
        cand_b[nc] = b;
        cand_t[nc++] = t;
    }
    for (int k = 0; k < 2; k++) {
        double t = k == 0 ? 0.0 : edge_len;
        double bnum = 0.0;
        for (int i = 0; i < n; i++) bnum += w[i] * (y[i] - s[i] * t);
        double b = sw ? bnum / sw : 0.0;
        if (b < 0.0) b = 0.0;
        cand_b[nc] = b;
        cand_t[nc++] = t;
        cand_b[nc] = 0.0;
        cand_t[nc++] = t;
    }
    if (sss > 0.0) {
        double t = ssy / sss;
        if (t < 0.0) t = 0.0;
        if (t > edge_len) t = edge_len;
        cand_b[nc] = 0.0;
        cand_t[nc++] = t;
    }
    for (int c = 0; c < nc; c++) {
        double rss = 0.0;
        for (int i = 0; i < n; i++) {
            double r = y[i] - cand_b[c] - s[i] * cand_t[c];
            rss += w[i] * r * r;
        }
        double mean = sw ? rss / sw : DBL_MAX;
        if (mean < best_rss) {
            best_rss = mean;
            best_b = cand_b[c];
            best_t = cand_t[c];
        }
    }
    *rss_out = best_rss;
    *pendant_out = best_b;
    *distal_child_out = best_t;
}

static int cmp_wls(const void *a, const void *b) {
    const placement_t *x = (const placement_t *)a;
    const placement_t *y = (const placement_t *)b;
    if (x->rss < y->rss) return -1;
    if (x->rss > y->rss) return 1;
    if (x->clade_size != y->clade_size) return x->clade_size - y->clade_size;
    return x->edge_child - y->edge_child;
}

static int cmp_sparse(const void *a, const void *b) {
    const placement_t *x = (const placement_t *)a;
    const placement_t *y = (const placement_t *)b;
    int x_ok = x->sparse_informative_ctx > 0;
    int y_ok = y->sparse_informative_ctx > 0;
    if (x_ok != y_ok) return y_ok - x_ok;
    if (x->sparse_score_per_ctx < y->sparse_score_per_ctx) return 1;
    if (x->sparse_score_per_ctx > y->sparse_score_per_ctx) return -1;
    if (x->sparse_score < y->sparse_score) return 1;
    if (x->sparse_score > y->sparse_score) return -1;
    return x->edge_child - y->edge_child;
}

static int cmp_final(const void *a, const void *b) {
    const placement_t *x = (const placement_t *)a;
    const placement_t *y = (const placement_t *)b;
    if (x->final_key < y->final_key) return -1;
    if (x->final_key > y->final_key) return 1;
    if (x->wls_rank != y->wls_rank) return x->wls_rank - y->wls_rank;
    if (x->sparse_rank != y->sparse_rank) return x->sparse_rank - y->sparse_rank;
    return x->edge_child - y->edge_child;
}

static void active_mask_for_query(const query_t *q, int n_refs, uint64_t *mask, int n_words) {
    memset(mask, 0, (size_t)n_words * sizeof(uint64_t));
    for (int r = 0; r < n_refs; r++) {
        if (q->has_dist[r]) bit_set(mask, r);
    }
}

static int nearest_ref_for_query(const query_t *q, int n_refs, double *dist_out) {
    int best = -1;
    double best_d = DBL_MAX;
    for (int r = 0; r < n_refs; r++) {
        if (!q->has_dist[r]) continue;
        if (q->dist[r] < best_d) {
            best_d = q->dist[r];
            best = r;
        }
    }
    if (dist_out) *dist_out = best_d;
    return best;
}

static void top_near_mask(const query_t *q, int n_refs, int k, uint64_t *mask, int n_words) {
    memset(mask, 0, (size_t)n_words * sizeof(uint64_t));
    if (k <= 0) return;
    unsigned char *used = (unsigned char *)xcalloc((size_t)n_refs, 1);
    for (int rank = 0; rank < k; rank++) {
        int best = -1;
        double best_d = DBL_MAX;
        for (int r = 0; r < n_refs; r++) {
            if (!q->has_dist[r] || used[r]) continue;
            if (q->dist[r] < best_d) {
                best_d = q->dist[r];
                best = r;
            }
        }
        if (best < 0) break;
        used[best] = 1;
        bit_set(mask, best);
    }
    free(used);
}

static int top_near_list(const query_t *q, int n_refs, int k, int *ids) {
    if (k <= 0) return 0;
    unsigned char *used = (unsigned char *)xcalloc((size_t)n_refs, 1);
    int n = 0;
    for (int rank = 0; rank < k; rank++) {
        int best = -1;
        double best_d = DBL_MAX;
        for (int r = 0; r < n_refs; r++) {
            if (!q->has_dist[r] || used[r]) continue;
            if (q->dist[r] < best_d) {
                best_d = q->dist[r];
                best = r;
            }
        }
        if (best < 0) break;
        used[best] = 1;
        ids[n++] = best;
    }
    free(used);
    return n;
}

static size_t lower_bound_ctx(const sorted_rec_t *arr, size_t n, uint64_t key) {
    size_t lo = 0;
    size_t hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (arr[mid].ctxgid < key) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static void compute_sparse_for_query(place_ctx_t *ctx, const query_t *q,
                                     placement_t *placements, int n_place) {
    if (!ctx->sparse_enabled || n_place <= 0) return;
    int sketch_index = q->sketch_index;
    if (sketch_index < 0 || (size_t)(sketch_index + 1) >= ctx->n_qidx) {
        die("no sparse sketch index is bound to query sample: %s", q->name);
    }
    int n_refs = ctx->refs->n;
    uint64_t start = ctx->qidx[sketch_index];
    uint64_t end = ctx->qidx[sketch_index + 1];
    if (end > ctx->n_qcomb) end = ctx->n_qcomb;

    int nearest_limit = ctx->sparse_nearest < n_refs ? ctx->sparse_nearest : n_refs;
    int *nearest_ids = (int *)xcalloc((size_t)nearest_limit, sizeof(int));
    int n_near = top_near_list(q, n_refs, nearest_limit, nearest_ids);
    int local_words = (n_near + 63) / 64;
    if (local_words <= 0) local_words = 1;
    size_t local_word_bytes = (size_t)local_words * sizeof(uint64_t);
    uint64_t local_ref_words[STACK_WORDS] = {0};
    uint64_t local_match_words[STACK_WORDS] = {0};
    int use_stack_words = local_words <= STACK_WORDS;
    uint64_t *ref_words = use_stack_words ? local_ref_words : (uint64_t *)xcalloc((size_t)local_words, sizeof(uint64_t));
    uint64_t *match_words = use_stack_words ? local_match_words : (uint64_t *)xcalloc((size_t)local_words, sizeof(uint64_t));
    uint64_t *edge_local_words = (uint64_t *)xcalloc((size_t)n_place * (size_t)local_words, sizeof(uint64_t));
    int *gid_to_local = (int *)xmalloc((size_t)n_refs * sizeof(int));
    for (int r = 0; r < n_refs; r++) gid_to_local[r] = -1;
    for (int i = 0; i < n_near; i++) gid_to_local[nearest_ids[i]] = i;

    for (int p = 0; p < n_place; p++) {
        uint64_t *local_edge = edge_local_words + (size_t)p * (size_t)local_words;
        uint64_t *edge_bits = ctx->node_bits + (size_t)placements[p].edge_child * (size_t)ctx->n_words;
        for (int i = 0; i < n_near; i++) {
            int gid = nearest_ids[i];
            int in_child_clade = bit_is_set(edge_bits, gid);
            int in_scoring_side = placements[p].sparse_uses_complement ? !in_child_clade : in_child_clade;
            if (in_scoring_side) bit_set(local_edge, i);
        }
    }

    int used_ctx = 0;
    for (uint64_t pos = start; pos < end; pos++) {
        uint64_t key = ctx->qcomb[pos];
        uint64_t ctx_id = key >> OBJ_BITS;
        uint32_t qobj = (uint32_t)(key & OBJ_MASK);
        size_t lo = lower_bound_ctx(ctx->sorted, ctx->n_sorted, ctx_id << GID_BITS);
        size_t hi = lower_bound_ctx(ctx->sorted, ctx->n_sorted, (ctx_id + 1) << GID_BITS);
        if (hi <= lo) continue;

        memset(ref_words, 0, local_word_bytes);
        memset(match_words, 0, local_word_bytes);
        for (size_t i = lo; i < hi; i++) {
            int gid = (int)(ctx->sorted[i].ctxgid & GID_MASK);
            if (gid < 0 || gid >= n_refs) continue;
            int local_id = gid_to_local[gid];
            if (local_id < 0) continue;
            bit_set(ref_words, local_id);
            if (ctx->sorted[i].obj == qobj) bit_set(match_words, local_id);
        }
        int n_total = popcount_words(ref_words, local_words);
        if (n_total < 2) continue;
        int m_total = popcount_words(match_words, local_words);
        if (m_total == 0 || m_total == n_total) continue;
        used_ctx++;

        for (int p = 0; p < n_place; p++) {
            uint64_t *local_edge = edge_local_words + (size_t)p * (size_t)local_words;
            int n_in = popcount_and(local_edge, ref_words, local_words);
            int n_out = n_total - n_in;
            if (n_in <= 0 || n_out <= 0) continue;
            int m_in = popcount_and(local_edge, match_words, local_words);
            int m_out = m_total - m_in;
            double p_in = ((double)m_in + 0.5) / ((double)n_in + 1.0);
            double p_out = ((double)m_out + 0.5) / ((double)n_out + 1.0);
            placements[p].sparse_score += log(p_in) - log(p_out);
            placements[p].sparse_informative_ctx++;
            placements[p].sparse_observed_ctx++;
        }
    }
    for (int p = 0; p < n_place; p++) {
        placements[p].query_ctx = end - start;
        placements[p].sparse_used_ctx = used_ctx;
        placements[p].sparse_score_per_ctx = placements[p].sparse_informative_ctx > 0
            ? placements[p].sparse_score / (double)placements[p].sparse_informative_ctx
            : -INFINITY;
    }

    if (!use_stack_words) {
        free(ref_words);
        free(match_words);
    }
    free(edge_local_words);
    free(gid_to_local);
    free(nearest_ids);
}

static char *node_label(const tree_t *tree, int node, char *buf, size_t n) {
    if (tree->nodes[node].name && tree->nodes[node].name[0]) return tree->nodes[node].name;
    snprintf(buf, n, "internal_%d", node);
    return buf;
}

static int newick_label_needs_quote(const char *s) {
    if (!s || !s[0]) return 0;
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
        if (isspace(*p) || strchr("'():,;[]{}", *p)) return 1;
    }
    return 0;
}

static void append_newick_label(strbuf_t *b, const char *s) {
    if (!s || !s[0]) return;
    if (!newick_label_needs_quote(s)) {
        strbuf_appendf(b, "%s", s);
        return;
    }
    strbuf_appendf(b, "'");
    for (const char *p = s; *p; p++) {
        if (*p == '\'') strbuf_appendf(b, "''");
        else strbuf_appendf(b, "%c", *p);
    }
    strbuf_appendf(b, "'");
}

static void append_jplace_tree_rec(strbuf_t *b, const tree_t *tree, int node) {
    int child = tree->nodes[node].first_child;
    if (child >= 0) {
        strbuf_appendf(b, "(");
        int first = 1;
        for (; child >= 0; child = tree->nodes[child].next_sibling) {
            if (!first) strbuf_appendf(b, ",");
            append_jplace_tree_rec(b, tree, child);
            first = 0;
        }
        strbuf_appendf(b, ")");
    }
    append_newick_label(b, tree->nodes[node].name);
    if (node != tree->root) {
        strbuf_appendf(b, ":%.10g{%d}", nonnegative_branch(tree->nodes[node].branch), node);
    }
}

static void append_json_string(strbuf_t *b, const char *s) {
    strbuf_appendf(b, "\"");
    if (s) {
        for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
            unsigned char c = *p;
            if (c == '\\' || c == '\"') strbuf_appendf(b, "\\%c", c);
            else if (c == '\b') strbuf_appendf(b, "\\b");
            else if (c == '\f') strbuf_appendf(b, "\\f");
            else if (c == '\n') strbuf_appendf(b, "\\n");
            else if (c == '\r') strbuf_appendf(b, "\\r");
            else if (c == '\t') strbuf_appendf(b, "\\t");
            else if (c < 0x20) strbuf_appendf(b, "\\u%04x", c);
            else strbuf_appendf(b, "%c", c);
        }
    }
    strbuf_appendf(b, "\"");
}

static char *build_invocation_string(int argc, char **argv) {
    strbuf_t b = {0};
    for (int i = 0; i < argc; i++) {
        if (i) strbuf_appendf(&b, " ");
        strbuf_appendf(&b, "%s", argv[i] ? argv[i] : "");
    }
    return b.data ? b.data : xstrdup("");
}

static double finite_or_zero(double x) {
    return isfinite(x) ? x : 0.0;
}

static void write_jplace_output(const char *path, const tree_t *tree, const queries_t *queries,
                                placement_t **tops, const int *n_tops, int argc, char **argv,
                                int negative_branches) {
    strbuf_t treebuf = {0};
    append_jplace_tree_rec(&treebuf, tree, tree->root);
    strbuf_appendf(&treebuf, ";");
    char *invocation = build_invocation_string(argc, argv);

    FILE *out = fopen(path, "w");
    if (!out) die("open jplace output failed for %s: %s", path, strerror(errno));
    fprintf(out, "{\n");
    fprintf(out, "  \"fields\": [\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\", \"pendant_length\"],\n");
    fprintf(out, "  \"metadata\": {\n");
    fprintf(out, "    \"invocation\": ");
    strbuf_t inv_json = {0};
    append_json_string(&inv_json, invocation);
    fputs(inv_json.data ? inv_json.data : "\"\"", out);
    fprintf(out, ",\n    \"program\": \"kssd3a place\",\n");
    fprintf(out, "    \"likelihood_note\": \"negative weighted least-squares residual, not a phylogenetic likelihood\",\n");
    fprintf(out, "    \"like_weight_ratio_note\": \"relative WLS support among reported candidate placements\",\n");
    fprintf(out, "    \"branch_lengths\": \"nonnegative p_dist lengths used by kssd3a place\",\n");
    fprintf(out, "    \"negative_branches_clamped\": %d\n", negative_branches);
    fprintf(out, "  },\n");
    fprintf(out, "  \"placements\": [\n");
    for (int qi = 0; qi < queries->n; qi++) {
        if (qi) fprintf(out, ",\n");
        fprintf(out, "    {\"n\": [");
        strbuf_t q_json = {0};
        append_json_string(&q_json, queries->items[qi].name);
        fputs(q_json.data ? q_json.data : "\"\"", out);
        free(q_json.data);
        fprintf(out, "], \"p\": [");
        double lwr_sum = 0.0;
        for (int k = 0; k < n_tops[qi]; k++) {
            lwr_sum += finite_or_zero(tops[qi][k].wls_rel_weight);
        }
        for (int k = 0; k < n_tops[qi]; k++) {
            placement_t *p = &tops[qi][k];
            double lwr = lwr_sum > 0.0 ? finite_or_zero(p->wls_rel_weight) / lwr_sum : (k == 0 ? 1.0 : 0.0);
            if (k) fprintf(out, ", ");
            fprintf(out, "[%d, %.12g, %.12g, %.12g, %.12g]",
                    p->edge_num,
                    -finite_or_zero(p->rss),
                    lwr,
                    finite_or_zero(p->distal_from_child),
                    finite_or_zero(p->pendant_length));
        }
        fprintf(out, "]}");
    }
    fprintf(out, "\n  ],\n  \"tree\": ");
    strbuf_t tree_json = {0};
    append_json_string(&tree_json, treebuf.data ? treebuf.data : ";");
    fputs(tree_json.data ? tree_json.data : "\";\"", out);
    fprintf(out, ",\n  \"version\": 3\n}\n");
    fclose(out);
    free(inv_json.data);
    free(tree_json.data);
    free(treebuf.data);
    free(invocation);
}

static void fill_depths_rec(const tree_t *tree, int node, int depth, double dist,
                            int *depths, double *root_dist) {
    depths[node] = depth;
    root_dist[node] = dist;
    for (int child = tree->nodes[node].first_child; child >= 0; child = tree->nodes[child].next_sibling) {
        fill_depths_rec(tree, child, depth + 1,
                        dist + nonnegative_branch(tree->nodes[child].branch),
                        depths, root_dist);
    }
}

static int lca_node(const tree_t *tree, const int *depths, int a, int b) {
    while (a >= 0 && b >= 0 && depths[a] > depths[b]) a = tree->nodes[a].parent;
    while (a >= 0 && b >= 0 && depths[b] > depths[a]) b = tree->nodes[b].parent;
    while (a >= 0 && b >= 0 && a != b) {
        a = tree->nodes[a].parent;
        b = tree->nodes[b].parent;
    }
    return a == b ? a : -1;
}

static double node_path_distance(const tree_t *tree, const int *depths, const double *root_dist,
                                 int a, int b) {
    int lca = lca_node(tree, depths, a, b);
    if (lca < 0) return NAN;
    return root_dist[a] + root_dist[b] - 2.0 * root_dist[lca];
}

static double placement_backbone_distance(const tree_t *tree, const int *depths, const double *root_dist,
                                          const placement_t *a, const placement_t *b) {
    if (!a || !b || a->edge_child < 0 || b->edge_child < 0) return NAN;
    int ac = a->edge_child;
    int bc = b->edge_child;
    int ap = tree->nodes[ac].parent;
    int bp = tree->nodes[bc].parent;
    if (ap < 0 || bp < 0) return NAN;
    if (ac == bc) return fabs(a->distal_from_parent - b->distal_from_parent);

    double candidates[4];
    candidates[0] = a->distal_from_parent + node_path_distance(tree, depths, root_dist, ap, bp) + b->distal_from_parent;
    candidates[1] = a->distal_from_parent + node_path_distance(tree, depths, root_dist, ap, bc) + b->distal_from_child;
    candidates[2] = a->distal_from_child + node_path_distance(tree, depths, root_dist, ac, bp) + b->distal_from_parent;
    candidates[3] = a->distal_from_child + node_path_distance(tree, depths, root_dist, ac, bc) + b->distal_from_child;
    double best = DBL_MAX;
    for (int i = 0; i < 4; i++) {
        if (isfinite(candidates[i]) && candidates[i] < best) best = candidates[i];
    }
    return best == DBL_MAX ? NAN : best;
}

static void write_pairwise_output(const char *path, const tree_t *tree, const queries_t *queries,
                                  placement_t **tops, const int *n_tops, char **warnings) {
    int *depths = (int *)xcalloc((size_t)tree->n, sizeof(int));
    double *root_dist = (double *)xcalloc((size_t)tree->n, sizeof(double));
    fill_depths_rec(tree, tree->root, 0, 0.0, depths, root_dist);

    FILE *out = fopen(path, "w");
    if (!out) die("open pairwise output failed for %s: %s", path, strerror(errno));
    fprintf(out, "query1\tquery2\tstatus\tplacement_p_dist\tbackbone_path_p_dist\tpendant_length1\tpendant_length2\tedge_num1\tedge_num2\tedge_child1\tedge_child2\tdistal_from_parent1\tdistal_from_parent2\tdistal_from_child1\tdistal_from_child2\tedge_length1\tedge_length2\twarnings1\twarnings2\n");
    for (int i = 0; i < queries->n; i++) {
        for (int j = i + 1; j < queries->n; j++) {
            const char *w1 = warnings && warnings[i] ? warnings[i] : "";
            const char *w2 = warnings && warnings[j] ? warnings[j] : "";
            if (n_tops[i] <= 0 || n_tops[j] <= 0) {
                fprintf(out, "%s\t%s\tmissing_placement\tNA\tNA\tNA\tNA\t-1\t-1\t-1\t-1\tNA\tNA\tNA\tNA\tNA\tNA\t%s\t%s\n",
                        queries->items[i].name, queries->items[j].name, w1, w2);
                continue;
            }
            placement_t *a = &tops[i][0];
            placement_t *b = &tops[j][0];
            double backbone = placement_backbone_distance(tree, depths, root_dist, a, b);
            double total = backbone + a->pendant_length + b->pendant_length;
            fprintf(out, "%s\t%s\tok\t%.12g\t%.12g\t%.12g\t%.12g\t%d\t%d\t%d\t%d\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%s\t%s\n",
                    queries->items[i].name,
                    queries->items[j].name,
                    total,
                    backbone,
                    a->pendant_length,
                    b->pendant_length,
                    a->edge_num,
                    b->edge_num,
                    a->edge_child,
                    b->edge_child,
                    a->distal_from_parent,
                    b->distal_from_parent,
                    a->distal_from_child,
                    b->distal_from_child,
                    a->edge_length,
                    b->edge_length,
                    w1,
                    w2);
        }
    }
    fclose(out);
    free(depths);
    free(root_dist);
}

static const char *rank_mode_name(rank_mode_t mode) {
    if (mode == RANK_SPARSE) return "sparse";
    if (mode == RANK_SUM) return "rank_sum";
    return "wls";
}

static char *place_one_query(place_ctx_t *ctx, int qi) {
    const query_t *q = &ctx->queries->items[qi];
    int n_refs = ctx->refs->n;
    int n_words = ctx->n_words;
    uint64_t *active = (uint64_t *)xcalloc((size_t)n_words, sizeof(uint64_t));
    uint64_t *near = (uint64_t *)xcalloc((size_t)n_words, sizeof(uint64_t));
    active_mask_for_query(q, n_refs, active, n_words);
    int active_count = popcount_words(active, n_words);
    if (ctx->candidate_near > 0) top_near_mask(q, n_refs, ctx->candidate_near, near, n_words);

    double nearest_dist = 0.0;
    int nearest = nearest_ref_for_query(q, n_refs, &nearest_dist);
    placement_t *placements = (placement_t *)xcalloc((size_t)ctx->tree->n, sizeof(placement_t));
    int n_place = 0;
    double *y = (double *)xcalloc((size_t)n_refs, sizeof(double));
    double *s = (double *)xcalloc((size_t)n_refs, sizeof(double));
    double *w = (double *)xcalloc((size_t)n_refs, sizeof(double));

    for (int child = 0; child < ctx->tree->n; child++) {
        int parent = ctx->tree->nodes[child].parent;
        if (parent < 0) continue;
        uint64_t *bits = ctx->node_bits + (size_t)child * n_words;
        int clade_size = popcount_and(bits, active, n_words);
        if (clade_size <= 0 || clade_size >= active_count) continue;
        if (ctx->candidate_near > 0 && popcount_and(bits, near, n_words) <= 0) continue;
        double edge_len = nonnegative_branch(ctx->tree->nodes[child].branch);
        int n = 0;
        for (int r = 0; r < n_refs; r++) {
            if (!q->has_dist[r]) continue;
            if (bit_is_set(bits, r)) {
                y[n] = q->dist[r] - ctx->node_dist[(size_t)child * n_refs + r];
                s[n] = 1.0;
            } else {
                y[n] = q->dist[r] - (ctx->node_dist[(size_t)parent * n_refs + r] + edge_len);
                s[n] = -1.0;
            }
            w[n] = weight_for(q->dist[r], ctx->weight_mode);
            n++;
        }
        if (n < 3) continue;
        double rss = 0.0, pendant = 0.0, distal_child = 0.0;
        fit_edge(y, s, w, n, edge_len, &rss, &pendant, &distal_child);
        placements[n_place].edge_child = child;
        placements[n_place].edge_num = child;
        placements[n_place].clade_size = clade_size;
        placements[n_place].rss = rss;
        placements[n_place].pendant_length = pendant;
        placements[n_place].distal_from_child = distal_child;
        placements[n_place].distal_from_parent = edge_len - distal_child;
        placements[n_place].edge_length = edge_len;
        placements[n_place].sparse_rank = ctx->tree->n + 1;
        int full_clade_size = popcount_words(bits, n_words);
        placements[n_place].sparse_uses_complement = full_clade_size > n_refs / 2;
        placements[n_place].sparse_key_size = placements[n_place].sparse_uses_complement
            ? n_refs - full_clade_size
            : full_clade_size;
        n_place++;
    }

    qsort(placements, (size_t)n_place, sizeof(placement_t), cmp_wls);
    for (int i = 0; i < n_place; i++) placements[i].wls_rank = i + 1;
    compute_wls_confidence(ctx, placements, n_place);

    compute_sparse_for_query(ctx, q, placements, n_place);
    qsort(placements, (size_t)n_place, sizeof(placement_t), cmp_sparse);
    int sparse_rank = 1;
    for (int i = 0; i < n_place; i++) {
        if (placements[i].sparse_informative_ctx > 0) placements[i].sparse_rank = sparse_rank++;
        else placements[i].sparse_rank = n_place + 1;
    }
    int wls_sparse_conflict = 0;
    if (ctx->sparse_enabled && n_place > 0) {
        int wls_top_idx = -1;
        int sparse_top_idx = -1;
        for (int i = 0; i < n_place; i++) {
            if (placements[i].wls_rank == 1) wls_top_idx = i;
            if (placements[i].sparse_rank == 1) sparse_top_idx = i;
        }
        if (wls_top_idx >= 0 && sparse_top_idx >= 0
            && placements[wls_top_idx].edge_child != placements[sparse_top_idx].edge_child
            && placements[wls_top_idx].sparse_rank > 10
            && placements[sparse_top_idx].wls_rank > 10) {
            wls_sparse_conflict = 1;
        }
    }
    for (int i = 0; i < n_place; i++) {
        if (ctx->rank_mode == RANK_SPARSE) placements[i].final_key = placements[i].sparse_rank;
        else if (ctx->rank_mode == RANK_SUM) placements[i].final_key = placements[i].wls_rank + placements[i].sparse_rank;
        else placements[i].final_key = placements[i].wls_rank;
    }
    qsort(placements, (size_t)n_place, sizeof(placement_t), cmp_final);

    strbuf_t out = {0};
    char warn[512] = "";
    if (active_count < 3) append_warning(warn, sizeof(warn), "too_few_ref_distances");
    if (n_place == 0) append_warning(warn, sizeof(warn), "no_candidate_edges");
    if (ctx->sparse_enabled && ctx->sparse_nearest < 2) append_warning(warn, sizeof(warn), "sparse_nearest_lt_2");
    if (n_place > 1 && placements[0].wls_top_weight > 0.0
        && placements[0].wls_second_weight / placements[0].wls_top_weight >= ctx->ambiguous_ratio) {
        append_warning(warn, sizeof(warn), "AMBIGUOUS_WLS");
    }
    if (ctx->sparse_enabled && n_place > 0
        && placements[0].sparse_used_ctx < ctx->min_sparse_informative) {
        append_warning(warn, sizeof(warn), "LOW_CONTEXT_SUPPORT");
    }
    if (wls_sparse_conflict) append_warning(warn, sizeof(warn), "WLS_SPARSE_CONFLICT");

    int limit = ctx->top_n < n_place ? ctx->top_n : n_place;
    const char *final_warn = warn[0] ? warn : (limit == 0 ? "no_placement" : "");
    ctx->n_top_placements[qi] = limit;
    ctx->warnings[qi] = xstrdup(final_warn);
    if (limit > 0) {
        ctx->top_placements[qi] = (placement_t *)xmalloc((size_t)limit * sizeof(placement_t));
        memcpy(ctx->top_placements[qi], placements, (size_t)limit * sizeof(placement_t));
    }
    for (int rank = 0; rank < limit; rank++) {
        placement_t *p = &placements[rank];
        char edge_buf[64];
        const char *edge_name = node_label(ctx->tree, p->edge_child, edge_buf, sizeof(edge_buf));
        strbuf_appendf(
            &out,
            "%s\t%d\t%s\t%d\t%s\t%d\t%s\t%.8g\t%d\t%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%d\t%d\t%d\t%.8g\t%.8g\t%d\t%d\t%" PRIu64 "\t%d\t%s\t%s\n",
            q->name,
            rank + 1,
            rank_mode_name(ctx->rank_mode),
            p->edge_child,
            edge_name,
            p->clade_size,
            nearest >= 0 ? ctx->refs->samples[nearest] : "",
            nearest >= 0 ? nearest_dist : NAN,
            n_place,
            p->wls_rank,
            p->rss,
            p->wls_delta,
            p->wls_delta_ratio,
            p->wls_rel_weight,
            p->wls_top_weight,
            p->wls_second_weight,
            p->wls_entropy,
            p->distal_from_parent,
            p->distal_from_child,
            p->pendant_length,
            p->edge_length,
            p->sparse_rank,
            p->sparse_key_size,
            p->sparse_uses_complement,
            p->sparse_score,
            p->sparse_score_per_ctx,
            p->sparse_informative_ctx,
            p->sparse_observed_ctx,
            p->query_ctx,
            p->sparse_used_ctx,
            ctx->weight_mode,
            final_warn
        );
    }
    if (limit == 0) {
        strbuf_appendf(&out, "%s\t0\t%s\t-1\t\t0\t\tNA\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\tNA\tNA\t0\t0\t0\t0\t%s\t%s\n",
                      q->name, rank_mode_name(ctx->rank_mode), ctx->weight_mode,
                      final_warn);
    }
    free(active);
    free(near);
    free(placements);
    free(y);
    free(s);
    free(w);
    return out.data ? out.data : xstrdup("");
}

static void *worker_main(void *argp) {
    worker_arg_t *arg = (worker_arg_t *)argp;
    place_ctx_t *ctx = arg->ctx;
    for (int qi = arg->tid; qi < ctx->queries->n; qi += ctx->n_threads) {
        ctx->blocks[qi] = place_one_query(ctx, qi);
    }
    return NULL;
}

static rank_mode_t parse_rank_mode(const char *s) {
    if (strcmp(s, "wls") == 0) return RANK_WLS;
    if (strcmp(s, "sparse") == 0) return RANK_SPARSE;
    if (strcmp(s, "rank_sum") == 0) return RANK_SUM;
    die("unknown --rank-by: %s", s);
    return RANK_WLS;
}

static void usage(FILE *fp) {
    fprintf(fp,
        "Usage: kssd3a place --tree TREE.nwk --idmap ref.idmap.tsv "
        "--distances DISTANCES --out placement.tsv [options]\n\n"
        "Options:\n"
        "  --ref-sketch DIR       reference -T sketch with sortedcomb_ctxgid64obj32\n"
        "  --query-sketch DIR     query -T sketch with comblco.index and comblco\n"
        "  --query-sketch-list FILE\n"
        "                         one query -T sketch directory per line; mutually exclusive\n"
        "                         with --query-sketch and combined in memory\n"
        "  --query-paths FILE     optional legacy query-label assertion; labels must\n"
        "                         exactly match query sketch metadata order\n"
        "  --distance-format FMT  auto|ani|matrix|phylip [auto]\n"
        "                         ani expects Qry/Ref/Distance columns; matrix expects\n"
        "                         tabular kssd3a matrix full output; phylip expects a\n"
        "                         square PHYLIP matrix with query and reference labels\n"
        "  --weight-mode MODE     equal|inv_d|inv_d2|inv_d_sqrt|inv_d_cap|close_exp [inv_d_cap]\n"
        "  --rank-by MODE         wls|sparse|rank_sum [wls]\n"
        "  --top N                output top N candidate edges per query [5]\n"
        "  --candidate-near N     WLS: score only edges containing at least one top-N nearest ref [0=all]\n"
        "  --sparse-nearest N     sparse object scoring uses top-N nearest refs [200]\n"
        "  --confidence-beta X    relative WLS weight temperature on normalized RSS gaps [0.5]\n"
        "  --ambiguous-ratio X    warn if second WLS relative weight / top weight >= X [0.8]\n"
        "  --min-sparse-informative N warn if sparse-used contexts below N [20]\n"
        "  --jplace-out FILE      write jplace output with edge-numbered tree\n"
        "  --pairwise-out FILE    write rank-1 query placement p-distance table\n"
        "  -o, --out FILE         placement TSV output path\n  -p, --threads N        query-level threads [1]\n"
        "  --help                 show this help\n");
}

int cmd_place(struct argp_state *state) {
    int argc = state->argc - state->next + 1;
    char **argv = &state->argv[state->next - 1];
    const char *tree_path = NULL;
    const char *idmap_path = NULL;
    const char *dist_path = NULL;
    const char *out_path = NULL;
    const char *jplace_out = NULL;
    const char *pairwise_out = NULL;
    const char *ref_sketch = NULL;
    const char *query_sketch = NULL;
    const char *query_sketch_list = NULL;
    const char *query_paths = NULL;
    distance_format_t distance_format = DISTANCE_FORMAT_AUTO;
    const char *weight_mode = "inv_d_cap";
    rank_mode_t rank_mode = RANK_WLS;
    int top_n = 5;
    int threads = 1;
    int candidate_near = 0;
    int sparse_nearest = 200;
    double confidence_beta = 0.5;
    double ambiguous_ratio = 0.8;
    int min_sparse_informative = 20;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--tree") == 0 && i + 1 < argc) tree_path = argv[++i];
        else if (strcmp(argv[i], "--idmap") == 0 && i + 1 < argc) idmap_path = argv[++i];
        else if (strcmp(argv[i], "--distances") == 0 && i + 1 < argc) dist_path = argv[++i];
        else if ((strcmp(argv[i], "--out") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) out_path = argv[++i];
        else if (strcmp(argv[i], "--jplace-out") == 0 && i + 1 < argc) jplace_out = argv[++i];
        else if (strcmp(argv[i], "--pairwise-out") == 0 && i + 1 < argc) pairwise_out = argv[++i];
        else if (strcmp(argv[i], "--ref-sketch") == 0 && i + 1 < argc) ref_sketch = argv[++i];
        else if (strcmp(argv[i], "--query-sketch") == 0 && i + 1 < argc) query_sketch = argv[++i];
        else if (strcmp(argv[i], "--query-sketch-list") == 0 && i + 1 < argc) query_sketch_list = argv[++i];
        else if (strcmp(argv[i], "--query-paths") == 0 && i + 1 < argc) query_paths = argv[++i];
        else if (strcmp(argv[i], "--distance-format") == 0 && i + 1 < argc) distance_format = parse_distance_format(argv[++i]);
        else if (strcmp(argv[i], "--weight-mode") == 0 && i + 1 < argc) weight_mode = argv[++i];
        else if (strcmp(argv[i], "--rank-by") == 0 && i + 1 < argc) rank_mode = parse_rank_mode(argv[++i]);
        else if (strcmp(argv[i], "--top") == 0 && i + 1 < argc) top_n = atoi(argv[++i]);
        else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-p") == 0) && i + 1 < argc) threads = atoi(argv[++i]);
        else if (strcmp(argv[i], "--candidate-near") == 0 && i + 1 < argc) candidate_near = atoi(argv[++i]);
        else if (strcmp(argv[i], "--sparse-nearest") == 0 && i + 1 < argc) sparse_nearest = atoi(argv[++i]);
        else if (strcmp(argv[i], "--confidence-beta") == 0 && i + 1 < argc) confidence_beta = atof(argv[++i]);
        else if (strcmp(argv[i], "--ambiguous-ratio") == 0 && i + 1 < argc) ambiguous_ratio = atof(argv[++i]);
        else if (strcmp(argv[i], "--min-sparse-informative") == 0 && i + 1 < argc) min_sparse_informative = atoi(argv[++i]);
        else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-?") == 0) {
            usage(stdout);
            state->next += argc - 1;
            return 0;
        } else {
            usage(stderr);
            return 2;
        }
    }
    if (!tree_path || !idmap_path || !dist_path || !out_path) {
        usage(stderr);
        return 2;
    }
    if (query_sketch && query_sketch_list) {
        die("--query-sketch and --query-sketch-list are mutually exclusive");
    }
    int sparse_enabled = ref_sketch || query_sketch || query_sketch_list;
    if (sparse_enabled && (!ref_sketch || (!query_sketch && !query_sketch_list))) {
        die("--ref-sketch and one of --query-sketch/--query-sketch-list are required for sparse-object scoring");
    }
    if (query_paths && !sparse_enabled) {
        die("--query-paths requires --ref-sketch and a query sketch");
    }
    if (top_n <= 0) top_n = 1;
    if (threads <= 0) threads = 1;
    if (sparse_nearest <= 0) sparse_nearest = 1;
    if (confidence_beta < 0.0) confidence_beta = 0.0;
    if (ambiguous_ratio <= 0.0 || ambiguous_ratio > 1.0) ambiguous_ratio = 0.8;
    if (min_sparse_informative < 0) min_sparse_informative = 0;

    double t0 = now_seconds();
    tree_t tree = read_newick(tree_path);
    refs_t refs = read_idmap(idmap_path);
    map_tree_refs(&tree, &refs);
    int n_words = 0;
    uint64_t *node_bits = build_node_bits(&tree, refs.n, &n_words);
    double *node_dist = precompute_node_ref_dist(&tree, &refs);

    sorted_rec_t *sorted = NULL;
    uint64_t *qidx = NULL;
    uint64_t *qcomb = NULL;
    char **query_names = NULL;
    size_t n_query_names = 0;
    size_t n_sorted = 0, n_qidx = 0, n_qcomb = 0;
    if (sparse_enabled) {
        char *sorted_path = path_join(ref_sketch, "sortedcomb_ctxgid64obj32");
        sorted = (sorted_rec_t *)read_binary_file(sorted_path, sizeof(sorted_rec_t), &n_sorted);
        free(sorted_path);
        if (query_sketch_list) {
            load_query_sketch_list(query_sketch_list, &qidx, &n_qidx, &qcomb, &n_qcomb,
                                   &query_names, &n_query_names);
        } else {
            load_query_sketch_dir(query_sketch, &qidx, &n_qidx, &qcomb, &n_qcomb,
                                  &query_names, &n_query_names);
        }
        if (n_qidx != n_query_names + 1) {
            die("query sketch metadata/index sample count mismatch");
        }
    }
    queries_t queries = read_distances(dist_path, &refs, query_paths, query_names, n_query_names,
                                       distance_format);
    free_query_sketch_names(query_names, n_query_names);
    if (queries.n <= 0) die("no query rows loaded from %s", dist_path);
    if (sparse_enabled && n_qidx != (size_t)queries.n + 1) {
        die("query sketch index has %zu records; expected %d for query count", n_qidx, queries.n + 1);
    }
    if (sparse_enabled) {
        for (int i = 0; i < queries.n; i++) {
            if (queries.items[i].n_dist == 0) {
                die("distance input has no reference distances for query sketch sample: %s",
                    queries.items[i].name);
            }
        }
    }

    int neg_branches = count_negative_branches(&tree);
    double setup_sec = now_seconds() - t0;
    if (threads > queries.n) threads = queries.n;

    place_ctx_t ctx = {
        .tree = &tree,
        .refs = &refs,
        .queries = &queries,
        .n_words = n_words,
        .node_bits = node_bits,
        .node_dist = node_dist,
        .weight_mode = weight_mode,
        .rank_mode = rank_mode,
        .top_n = top_n,
        .candidate_near = candidate_near,
        .sparse_nearest = sparse_nearest,
        .n_threads = threads,
        .confidence_beta = confidence_beta,
        .ambiguous_ratio = ambiguous_ratio,
        .min_sparse_informative = min_sparse_informative,
        .sparse_enabled = sparse_enabled,
        .sorted = sorted,
        .n_sorted = n_sorted,
        .qidx = qidx,
        .n_qidx = n_qidx,
        .qcomb = qcomb,
        .n_qcomb = n_qcomb,
        .blocks = (char **)xcalloc((size_t)queries.n, sizeof(char *)),
        .top_placements = (placement_t **)xcalloc((size_t)queries.n, sizeof(placement_t *)),
        .n_top_placements = (int *)xcalloc((size_t)queries.n, sizeof(int)),
        .warnings = (char **)xcalloc((size_t)queries.n, sizeof(char *)),
    };

    double t1 = now_seconds();
    pthread_t *ths = (pthread_t *)xcalloc((size_t)threads, sizeof(pthread_t));
    worker_arg_t *args = (worker_arg_t *)xcalloc((size_t)threads, sizeof(worker_arg_t));
    for (int t = 0; t < threads; t++) {
        args[t].ctx = &ctx;
        args[t].tid = t;
        int rc = pthread_create(&ths[t], NULL, worker_main, &args[t]);
        if (rc != 0) die("pthread_create failed: %s", strerror(rc));
    }
    for (int t = 0; t < threads; t++) pthread_join(ths[t], NULL);
    double place_sec = now_seconds() - t1;

    FILE *out = fopen(out_path, "w");
    if (!out) die("open output failed for %s: %s", out_path, strerror(errno));
    fprintf(out, "#setup_seconds=%.6f\tplacement_seconds=%.6f\tqueries=%d\trefs=%d\tnodes=%d\tthreads=%d\tnegative_branches_clamped=%d\tsparse_enabled=%d\trank_by=%s\tsparse_nearest=%d\tcandidate_near=%d\tconfidence_beta=%.8g\tambiguous_ratio=%.8g\tmin_sparse_informative=%d\n",
            setup_sec, place_sec, queries.n, refs.n, tree.n, threads, neg_branches,
            sparse_enabled, rank_mode_name(rank_mode), sparse_nearest, candidate_near,
            confidence_beta, ambiguous_ratio, min_sparse_informative);
    fprintf(out, "query\trank\trank_by\tedge_child\tedge_label\tclade_size\tnearest_ref\tnearest_p_dist\tcandidate_count\twls_rank\twls_score\twls_delta\twls_delta_ratio\twls_rel_weight\twls_top_weight\twls_second_weight\twls_entropy\tdistal_from_parent\tdistal_from_child\tpendant_length\tedge_length\tsparse_rank\tsparse_key_size\tsparse_uses_complement\tsparse_score\tsparse_score_per_ctx\tsparse_informative_ctx\tsparse_observed_ctx\tquery_ctx\tsparse_used_ctx\tweight_mode\twarnings\n");
    for (int qi = 0; qi < queries.n; qi++) {
        if (ctx.blocks[qi]) fputs(ctx.blocks[qi], out);
    }
    fclose(out);

    if (jplace_out) {
        write_jplace_output(jplace_out, &tree, &queries, ctx.top_placements,
                            ctx.n_top_placements, argc, argv, neg_branches);
    }
    if (pairwise_out) {
        write_pairwise_output(pairwise_out, &tree, &queries, ctx.top_placements,
                              ctx.n_top_placements, ctx.warnings);
    }

    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    fprintf(stderr, "queries=%d refs=%d nodes=%d setup=%.3fs placement=%.3fs maxrss=%.1fMB negative_branches_clamped=%d sparse=%d rank_by=%s\n",
            queries.n, refs.n, tree.n, setup_sec, place_sec, ru.ru_maxrss / 1024.0,
            neg_branches, sparse_enabled, rank_mode_name(rank_mode));
    state->next += argc - 1;
    return 0;
}
