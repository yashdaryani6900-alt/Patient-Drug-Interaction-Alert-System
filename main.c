#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ------------------------------------------------------------------ */
/* Constants & types                                                  */
/* ------------------------------------------------------------------ */
#define NAME_LEN   40
#define NOTE_LEN   80
#define LINE_LEN   256
#define INIT_CAP   4
#define DB_FILE    "drugs.db"

typedef enum {
    SEV_NONE = 0,
    SEV_MINOR,
    SEV_MAJOR,
    SEV_CONTRA
} Severity;

typedef struct {
    int      otherId;           /* id of the interacting drug            */
    Severity severity;          /* minor / major / contraindicated       */
    char     note[NOTE_LEN];    /* short clinical note                   */
} Interaction;

typedef struct {
    int          id;
    char         name[NAME_LEN];
    char         drugClass[NAME_LEN];   /* therapeutic class             */
    float        minDose, maxDose;      /* dosage limits (mg/day)        */
    Interaction *interactions;          /* dynamic array                 */
    int          intCount, intCap;
} Drug;

typedef struct {
    Drug *items;                        /* dynamic array, sorted by id   */
    int   count, cap;
} DrugDB;

typedef struct {
    int   drugId;
    float dose;                         /* 0 = not given                 */
} PatientMed;

/* ------------------------------------------------------------------ */
/* Small utilities (strings & parsing)                                */
/* ------------------------------------------------------------------ */
static const char *severity_name(Severity s)
{
    switch (s) {
        case SEV_MINOR:  return "MINOR";
        case SEV_MAJOR:  return "MAJOR";
        case SEV_CONTRA: return "CONTRAINDICATED";
        default:         return "NONE";
    }
}

static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) s++;
    if (*s == '\0') return s;
    char *end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
    return s;
}

static int str_ieq(const char *a, const char *b)       /* case-insensitive == */
{
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return 0;
        a++; b++;
    }
    return *a == *b;
}

static void read_line(const char *prompt, char *buf, size_t n)
{
    printf("%s", prompt);
    fflush(stdout);
    if (!fgets(buf, (int)n, stdin)) { buf[0] = '\0'; return; }
    buf[strcspn(buf, "\r\n")] = '\0';
}

static int read_int(const char *prompt, int *out)
{
    char buf[LINE_LEN], *end;
    read_line(prompt, buf, sizeof buf);
    long v = strtol(buf, &end, 10);
    if (end == buf || *trim(end) != '\0') { printf("  Invalid number.\n"); return 0; }
    *out = (int)v;
    return 1;
}

static int read_float(const char *prompt, float *out)
{
    char buf[LINE_LEN], *end;
    read_line(prompt, buf, sizeof buf);
    float v = strtof(buf, &end);
    if (end == buf || *trim(end) != '\0') { printf("  Invalid number.\n"); return 0; }
    *out = v;
    return 1;
}

static void *xrealloc(void *p, size_t bytes)
{
    void *q = realloc(p, bytes);
    if (!q) { fprintf(stderr, "Out of memory\n"); exit(EXIT_FAILURE); }
    return q;
}

/* ------------------------------------------------------------------ */
/* Database: creation, search, maintenance (pointers + dynamic memory) */
/* ------------------------------------------------------------------ */
static void db_init(DrugDB *db)
{
    db->count = 0;
    db->cap   = INIT_CAP;
    db->items = xrealloc(NULL, (size_t)db->cap * sizeof(Drug));
}

static void db_free(DrugDB *db)
{
    for (int i = 0; i < db->count; i++) free(db->items[i].interactions);
    free(db->items);
    db->items = NULL;
    db->count = db->cap = 0;
}

/* Binary search by id (array is kept sorted). Returns index or -1. */
static int db_find_id(const DrugDB *db, int id)
{
    int lo = 0, hi = db->count - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        if (db->items[mid].id == id) return mid;
        if (db->items[mid].id < id)  lo = mid + 1;
        else                         hi = mid - 1;
    }
    return -1;
}

/* Linear search by name (case-insensitive). Returns index or -1. */
static int db_find_name(const DrugDB *db, const char *name)
{
    for (int i = 0; i < db->count; i++)
        if (str_ieq(db->items[i].name, name)) return i;
    return -1;
}

/* Add a drug, keeping the array sorted by id. Returns index or -1. */
static int db_add_drug(DrugDB *db, int id, const char *name, const char *cls,
                       float minDose, float maxDose)
{
    if (db_find_id(db, id) != -1 || db_find_name(db, name) != -1) {
        printf("  Drug id or name already exists.\n");
        return -1;
    }
    if (db->count == db->cap) {
        db->cap *= 2;
        db->items = xrealloc(db->items, (size_t)db->cap * sizeof(Drug));
    }
    int pos = 0;
    while (pos < db->count && db->items[pos].id < id) pos++;
    memmove(&db->items[pos + 1], &db->items[pos],
            (size_t)(db->count - pos) * sizeof(Drug));

    Drug *d = &db->items[pos];
    memset(d, 0, sizeof *d);
    d->id = id;
    strncpy(d->name, name, NAME_LEN - 1);
    strncpy(d->drugClass, cls, NAME_LEN - 1);
    d->minDose = minDose;
    d->maxDose = maxDose;
    d->intCap  = INIT_CAP;
    d->interactions = xrealloc(NULL, (size_t)d->intCap * sizeof(Interaction));
    db->count++;
    return pos;
}

/* Find an interaction inside one drug's list (linear search). */
static Interaction *drug_find_interaction(Drug *d, int otherId)
{
    for (int i = 0; i < d->intCount; i++)
        if (d->interactions[i].otherId == otherId) return &d->interactions[i];
    return NULL;
}

/* One-directional insert/update. */
static void drug_set_interaction(Drug *d, int otherId, Severity sev, const char *note)
{
    Interaction *it = drug_find_interaction(d, otherId);
    if (!it) {
        if (d->intCount == d->intCap) {
            d->intCap *= 2;
            d->interactions = xrealloc(d->interactions,
                                       (size_t)d->intCap * sizeof(Interaction));
        }
        it = &d->interactions[d->intCount++];
        it->otherId = otherId;
    }
    it->severity = sev;
    strncpy(it->note, note, NOTE_LEN - 1);
    it->note[NOTE_LEN - 1] = '\0';
}

static void drug_remove_interaction(Drug *d, int otherId)
{
    for (int i = 0; i < d->intCount; i++) {
        if (d->interactions[i].otherId == otherId) {
            d->interactions[i] = d->interactions[--d->intCount];
            return;
        }
    }
}

/* Interactions are symmetric: store on both drugs. */
static int db_link(DrugDB *db, int idA, int idB, Severity sev, const char *note)
{
    int a = db_find_id(db, idA), b = db_find_id(db, idB);
    if (a < 0 || b < 0 || a == b) { printf("  Invalid drug pair.\n"); return 0; }
    drug_set_interaction(&db->items[a], idB, sev, note);
    drug_set_interaction(&db->items[b], idA, sev, note);
    return 1;
}

static void db_unlink(DrugDB *db, int idA, int idB)
{
    int a = db_find_id(db, idA), b = db_find_id(db, idB);
    if (a >= 0) drug_remove_interaction(&db->items[a], idB);
    if (b >= 0) drug_remove_interaction(&db->items[b], idA);
}

static void db_delete_drug(DrugDB *db, int id)
{
    int pos = db_find_id(db, id);
    if (pos < 0) { printf("  Drug not found.\n"); return; }
    free(db->items[pos].interactions);
    memmove(&db->items[pos], &db->items[pos + 1],
            (size_t)(db->count - pos - 1) * sizeof(Drug));
    db->count--;
    for (int i = 0; i < db->count; i++)           /* clean dangling references */
        drug_remove_interaction(&db->items[i], id);
    printf("  Drug deleted.\n");
}

/* Highest severity recorded between two drugs (checks both lists). */
static Severity severity_between(const DrugDB *db, int idA, int idB, const char **noteOut)
{
    Severity best = SEV_NONE;
    const char *note = "";
    int a = db_find_id(db, idA), b = db_find_id(db, idB);
    if (a >= 0) {
        Interaction *it = drug_find_interaction(&db->items[a], idB);
        if (it && it->severity > best) { best = it->severity; note = it->note; }
    }
    if (b >= 0) {
        Interaction *it = drug_find_interaction(&db->items[b], idA);
        if (it && it->severity > best) { best = it->severity; note = it->note; }
    }
    if (noteOut) *noteOut = note;
    return best;
}

/* ------------------------------------------------------------------ */
/* Display                                                            */
/* ------------------------------------------------------------------ */
static void print_drug(const DrugDB *db, const Drug *d)
{
    printf("[%d] %-16s class: %-16s dose: %.1f-%.1f mg/day\n",
           d->id, d->name, d->drugClass, d->minDose, d->maxDose);
    for (int i = 0; i < d->intCount; i++) {
        int k = db_find_id(db, d->interactions[i].otherId);
        printf("      - %-16s %-16s %s\n",
               k >= 0 ? db->items[k].name : "(removed)",
               severity_name(d->interactions[i].severity),
               d->interactions[i].note);
    }
}

static void list_drugs(const DrugDB *db)
{
    if (db->count == 0) { printf("  Database is empty.\n"); return; }
    for (int i = 0; i < db->count; i++) print_drug(db, &db->items[i]);
}

/* ------------------------------------------------------------------ */
/* Patient medication parsing:  "warfarin:5, aspirin, 3:20"           */
/* token = name-or-id [: dose]                                         */
/* ------------------------------------------------------------------ */
static int parse_patient_meds(const DrugDB *db, char *line, PatientMed **out)
{
    int cap = INIT_CAP, n = 0;
    PatientMed *meds = xrealloc(NULL, (size_t)cap * sizeof(PatientMed));

    for (char *tok = strtok(line, ","); tok; tok = strtok(NULL, ",")) {
        char *name = trim(tok);
        float dose = 0.0f;
        char *colon = strchr(name, ':');
        if (colon) {
            *colon = '\0';
            dose = strtof(trim(colon + 1), NULL);
            name = trim(name);
        }
        if (*name == '\0') continue;

        char *end;
        long asId = strtol(name, &end, 10);
        int idx = (*end == '\0') ? db_find_id(db, (int)asId)
                                 : db_find_name(db, name);
        if (idx < 0) { printf("  Unknown drug '%s' - skipped.\n", name); continue; }

        int dup = 0;
        for (int i = 0; i < n; i++)
            if (meds[i].drugId == db->items[idx].id) dup = 1;
        if (dup) continue;

        if (n == cap) { cap *= 2; meds = xrealloc(meds, (size_t)cap * sizeof(PatientMed)); }
        meds[n].drugId = db->items[idx].id;
        meds[n].dose   = dose;
        n++;
    }
    *out = meds;
    return n;
}

/* ------------------------------------------------------------------ */
/* Core logic: interaction check, dose check, alternatives            */
/* ------------------------------------------------------------------ */
static void check_doses(const DrugDB *db, const PatientMed *meds, int n)
{
    for (int i = 0; i < n; i++) {
        if (meds[i].dose <= 0.0f) continue;
        const Drug *d = &db->items[db_find_id(db, meds[i].drugId)];
        if (meds[i].dose > d->maxDose)
            printf("  [DOSE] %s: %.1f mg exceeds max %.1f mg/day\n", d->name, meds[i].dose, d->maxDose);
        else if (meds[i].dose < d->minDose)
            printf("  [DOSE] %s: %.1f mg is below min %.1f mg/day\n", d->name, meds[i].dose, d->minDose);
    }
}

/* Pairwise scan; prints most severe first. Returns number of flagged pairs. */
static int scan_interactions(const DrugDB *db, const PatientMed *meds, int n,
                             int *needsReplacement)   /* needsReplacement[n] flags */
{
    int total = 0, counts[SEV_CONTRA + 1] = {0};

    for (int sev = SEV_CONTRA; sev >= SEV_MINOR; sev--) {
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                const char *note;
                Severity s = severity_between(db, meds[i].drugId, meds[j].drugId, &note);
                if ((int)s != sev) continue;

                const char *a = db->items[db_find_id(db, meds[i].drugId)].name;
                const char *b = db->items[db_find_id(db, meds[j].drugId)].name;
                printf("  [%-15s] %s + %s : %s\n", severity_name(s), a, b, note);
                counts[s]++; total++;
                if (s >= SEV_MAJOR) { needsReplacement[i] = 1; needsReplacement[j] = 1; }
            }
        }
    }
    printf("\n  Summary: %d contraindicated, %d major, %d minor (%d pair(s) scanned)\n",
           counts[SEV_CONTRA], counts[SEV_MAJOR], counts[SEV_MINOR], n * (n - 1) / 2);
    return total;
}

/* Same class, not already prescribed, and no interaction with the other meds. */
static int suggest_alternatives(const DrugDB *db, const PatientMed *meds, int n, int targetId)
{
    int t = db_find_id(db, targetId);
    if (t < 0) { printf("  Drug not found.\n"); return 0; }
    const Drug *target = &db->items[t];
    int found = 0;

    printf("  Alternatives for %s (class: %s):\n", target->name, target->drugClass);
    for (int i = 0; i < db->count; i++) {
        const Drug *c = &db->items[i];
        if (c->id == targetId || !str_ieq(c->drugClass, target->drugClass)) continue;

        int ok = 1;
        for (int k = 0; k < n && ok; k++) {
            if (meds[k].drugId == targetId) continue;      /* being replaced */
            if (meds[k].drugId == c->id) { ok = 0; break; } /* already taking */
            if (severity_between(db, c->id, meds[k].drugId, NULL) != SEV_NONE) ok = 0;
        }
        if (ok) {
            printf("    -> %s (id %d), dose %.1f-%.1f mg/day\n", c->name, c->id, c->minDose, c->maxDose);
            found++;
        }
    }
    if (!found) printf("    No interaction-free alternative in the formulary.\n");
    return found;
}

static void analyze_patient(const DrugDB *db)
{
    char line[LINE_LEN];
    printf("  Enter medications (name or id, optional ':dose', comma separated)\n");
    read_line("  e.g. warfarin:5, aspirin:75, 4 > ", line, sizeof line);

    PatientMed *meds;
    int n = parse_patient_meds(db, line, &meds);
    if (n < 2) {
        printf("  Need at least two valid medications to check interactions.\n");
        if (n == 1) check_doses(db, meds, n);
        free(meds);
        return;
    }

    int *flag = calloc((size_t)n, sizeof(int));
    if (!flag) { free(meds); return; }

    printf("\n=== INTERACTION REPORT ===\n");
    int found = scan_interactions(db, meds, n, flag);
    if (found == 0) printf("  No known interactions.\n");

    printf("\n=== DOSAGE CHECK ===\n");
    check_doses(db, meds, n);

    int any = 0;
    for (int i = 0; i < n; i++) any |= flag[i];
    if (any) {
        printf("\n=== SUGGESTED REPLACEMENTS (major / contraindicated only) ===\n");
        for (int i = 0; i < n; i++)
            if (flag[i]) suggest_alternatives(db, meds, n, meds[i].drugId);
    }
    free(flag);
    free(meds);
}

/* ------------------------------------------------------------------ */
/* Persistence (plain text, easy to extend)                            */
/* D|id|name|class|min|max|interactionCount                            */
/* I|otherId|severity|note                                             */
/* ------------------------------------------------------------------ */
static int db_save(const DrugDB *db, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) return 0;
    for (int i = 0; i < db->count; i++) {
        const Drug *d = &db->items[i];
        fprintf(f, "D|%d|%s|%s|%.2f|%.2f|%d\n", d->id, d->name, d->drugClass,
                d->minDose, d->maxDose, d->intCount);
        for (int k = 0; k < d->intCount; k++)
            fprintf(f, "I|%d|%d|%s\n", d->interactions[k].otherId,
                    (int)d->interactions[k].severity, d->interactions[k].note);
    }
    fclose(f);
    return 1;
}

static int db_load(DrugDB *db, const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char line[LINE_LEN], name[NAME_LEN], cls[NAME_LEN], note[NOTE_LEN];
    int id, cnt, other, sev, cur = -1;
    float mn, mx;

    while (fgets(line, sizeof line, f)) {
        line[strcspn(line, "\r\n")] = '\0';
        if (line[0] == 'D' &&
            sscanf(line + 2, "%d|%39[^|]|%39[^|]|%f|%f|%d", &id, name, cls, &mn, &mx, &cnt) == 6) {
            int idx = db_add_drug(db, id, name, cls, mn, mx);
            cur = (idx >= 0) ? id : -1;
        } else if (line[0] == 'I' && cur >= 0 &&
                   sscanf(line + 2, "%d|%d|%79[^\n]", &other, &sev, note) == 3) {
            drug_set_interaction(&db->items[db_find_id(db, cur)], other, (Severity)sev, note);
        }
    }
    fclose(f);
    return 1;
}

static void seed_database(DrugDB *db)
{
    db_add_drug(db, 1, "Warfarin",     "Anticoagulant", 1.0f, 10.0f);
    db_add_drug(db, 2, "Aspirin",      "Antiplatelet",  75.0f, 325.0f);
    db_add_drug(db, 3, "Clopidogrel",  "Antiplatelet",  75.0f, 75.0f);
    db_add_drug(db, 4, "Ibuprofen",    "NSAID",         200.0f, 1200.0f);
    db_add_drug(db, 5, "Paracetamol",  "Analgesic",     500.0f, 4000.0f);
    db_add_drug(db, 6, "Naproxen",     "NSAID",         250.0f, 1250.0f);
    db_add_drug(db, 7, "Sertraline",   "SSRI",          50.0f, 200.0f);
    db_add_drug(db, 8, "Tramadol",     "Opioid",        50.0f, 400.0f);
    db_add_drug(db, 9, "Phenelzine",   "MAOI",          15.0f, 90.0f);
    db_add_drug(db, 10, "Escitalopram","SSRI",          10.0f, 20.0f);
    db_add_drug(db, 11, "Rivaroxaban", "Anticoagulant", 10.0f, 20.0f);

    db_link(db, 1, 2,  SEV_MAJOR,  "Increased bleeding risk");
    db_link(db, 1, 4,  SEV_MAJOR,  "Increased bleeding risk, GI bleed");
    db_link(db, 1, 6,  SEV_MAJOR,  "Increased bleeding risk");
    db_link(db, 1, 5,  SEV_MINOR,  "High-dose paracetamol may raise INR");
    db_link(db, 2, 4,  SEV_MINOR,  "Ibuprofen may reduce aspirin effect");
    db_link(db, 2, 3,  SEV_MAJOR,  "Dual antiplatelet: bleeding risk");
    db_link(db, 7, 8,  SEV_MAJOR,  "Serotonin syndrome risk");
    db_link(db, 7, 9,  SEV_CONTRA, "Serotonin syndrome - never combine");
    db_link(db, 8, 9,  SEV_CONTRA, "Severe serotonin syndrome risk");
    db_link(db, 10, 9, SEV_CONTRA, "Serotonin syndrome - never combine");
    db_link(db, 7, 2,  SEV_MINOR,  "Slightly higher bleeding risk");
    db_link(db, 11, 1, SEV_CONTRA, "Duplicate anticoagulation");
    db_link(db, 11, 4, SEV_MAJOR,  "Increased bleeding risk");
}

/* ------------------------------------------------------------------ */
/* Menu actions                                                       */
/* ------------------------------------------------------------------ */
static int pick_drug(const DrugDB *db, const char *prompt)
{
    char buf[LINE_LEN];
    read_line(prompt, buf, sizeof buf);
    char *s = trim(buf), *end;
    long v = strtol(s, &end, 10);
    int idx = (*s && *end == '\0') ? db_find_id(db, (int)v) : db_find_name(db, s);
    if (idx < 0) { printf("  Drug not found.\n"); return -1; }
    return db->items[idx].id;
}

static void menu_add_drug(DrugDB *db)
{
    int id; float mn, mx; char name[LINE_LEN], cls[LINE_LEN];
    if (!read_int("  Drug id: ", &id)) return;
    read_line("  Name: ", name, sizeof name);
    read_line("  Therapeutic class: ", cls, sizeof cls);
    if (!read_float("  Min dose (mg/day): ", &mn)) return;
    if (!read_float("  Max dose (mg/day): ", &mx)) return;
    if (*trim(name) == '\0' || mn > mx) { printf("  Invalid data.\n"); return; }
    if (db_add_drug(db, id, trim(name), trim(cls), mn, mx) >= 0) printf("  Drug added.\n");
}

static void menu_update_dose(DrugDB *db)
{
    int id = pick_drug(db, "  Drug (id or name): ");
    if (id < 0) return;
    float mn, mx;
    if (!read_float("  New min dose: ", &mn) || !read_float("  New max dose: ", &mx)) return;
    if (mn > mx) { printf("  Min cannot exceed max.\n"); return; }
    Drug *d = &db->items[db_find_id(db, id)];
    d->minDose = mn; d->maxDose = mx;
    printf("  Dosage limits updated.\n");
}

static void menu_set_interaction(DrugDB *db)
{
    int a = pick_drug(db, "  First drug: ");   if (a < 0) return;
    int b = pick_drug(db, "  Second drug: ");  if (b < 0) return;
    int sev; char note[LINE_LEN];
    if (!read_int("  Severity (1=minor, 2=major, 3=contraindicated): ", &sev)) return;
    if (sev < SEV_MINOR || sev > SEV_CONTRA) { printf("  Invalid severity.\n"); return; }
    read_line("  Note: ", note, sizeof note);
    if (db_link(db, a, b, (Severity)sev, trim(note))) printf("  Interaction saved.\n");
}

static void menu_remove_interaction(DrugDB *db)
{
    int a = pick_drug(db, "  First drug: ");   if (a < 0) return;
    int b = pick_drug(db, "  Second drug: ");  if (b < 0) return;
    db_unlink(db, a, b);
    printf("  Interaction removed (if it existed).\n");
}

static void menu_alternatives(const DrugDB *db)
{
    int id = pick_drug(db, "  Drug to replace (id or name): ");
    if (id < 0) return;
    char line[LINE_LEN];
    read_line("  Patient's other meds (comma separated, blank for none): ", line, sizeof line);
    PatientMed *meds;
    int n = parse_patient_meds(db, line, &meds);
    suggest_alternatives(db, meds, n, id);
    free(meds);
}

static void print_menu(void)
{
    printf("\n===== DRUG INTERACTION ALERT SYSTEM =====\n"
           " 1. List drug database\n"
           " 2. Add drug\n"
           " 3. Update dosage limits\n"
           " 4. Add / update interaction\n"
           " 5. Remove interaction\n"
           " 6. Delete drug\n"
           " 7. Check patient medication list\n"
           " 8. Suggest alternatives for a drug\n"
           " 9. Save database\n"
           " 0. Exit\n");
}

/* ------------------------------------------------------------------ */
int main(void)
{
    DrugDB db;
    db_init(&db);
    if (db_load(&db, DB_FILE) && db.count > 0) printf("Loaded %d drugs from %s\n", db.count, DB_FILE);
    else { seed_database(&db); printf("Loaded built-in sample formulary (%d drugs).\n", db.count); }

    int choice = -1;
    while (choice != 0) {
        print_menu();
        if (!read_int("Choice: ", &choice)) { choice = -1; continue; }
        switch (choice) {
            case 1: list_drugs(&db); break;
            case 2: menu_add_drug(&db); break;
            case 3: menu_update_dose(&db); break;
            case 4: menu_set_interaction(&db); break;
            case 5: menu_remove_interaction(&db); break;
            case 6: { int id = pick_drug(&db, "  Drug to delete: "); if (id >= 0) db_delete_drug(&db, id); break; }
            case 7: analyze_patient(&db); break;
            case 8: menu_alternatives(&db); break;
            case 9: printf(db_save(&db, DB_FILE) ? "  Saved to %s\n" : "  Could not save %s\n", DB_FILE); break;
            case 0: break;
            default: printf("  Unknown option.\n");
        }
    }
    db_free(&db);
    printf("Goodbye.\n");
    return 0;
}
