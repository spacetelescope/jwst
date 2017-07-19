from docutils.nodes import figure, caption, Text, reference, raw, SkipNode, Element
from sphinx.roles import XRefRole


# Element classes

class page_ref(reference):
    pass

class num_ref(reference):
    pass


# Visit/depart functions

def skip_page_ref(self, node):
    raise SkipNode

def latex_visit_page_ref(self, node):
    self.body.append("\\pageref{%s:%s}" % (node['refdoc'], node['reftarget']))
    raise SkipNode

def latex_visit_num_ref(self, node):
    fields = node['reftarget'].split('#')
    if len(fields) > 1:
        label, target = fields
        ref_link = '%s:%s' % (node['refdoc'], target)
        latex = "\\hyperref[%s]{%s \\ref*{%s}}" % (ref_link, label, ref_link)
        self.body.append(latex)
    else:
        self.body.append('\\ref{%s:%s}' % (node['refdoc'], fields[0]))

    raise SkipNode


def doctree_read(app, doctree):
    # first generate figure numbers for each figure
    env = app.builder.env
    figid_docname_map = getattr(env, 'figid_docname_map', {})

    for figure_info in doctree.traverse(figure):
        for id in figure_info['ids']:
            figid_docname_map[id] = env.docname

    env.figid_docname_map = figid_docname_map


def doctree_resolved(app, doctree, docname):
    i = 1
    figids = {}
    for figure_info in doctree.traverse(figure):
        if app.builder.name != 'latex' and app.config.number_figures:
            for cap in figure_info.traverse(caption):
                cap[0] = Text("%s %d: %s" % (app.config.figure_caption_prefix, i, cap[0]))

        for id in figure_info['ids']:
            figids[id] = i

        i += 1


    # replace numfig nodes with links
    if app.builder.name != 'latex':
        for ref_info in doctree.traverse(num_ref):
            if '#' in ref_info['reftarget']:
                label, target = ref_info['reftarget'].split('#')
                labelfmt = label + " %d"
            else:
                labelfmt = '%d'
                target = ref_info['reftarget']

            if target not in figids:
                continue

            if app.builder.name == 'html':
                target_doc = app.builder.env.figid_docname_map[target]
                link = "%s#%s" % (app.builder.get_relative_uri(docname, target_doc),
                                  target)
                html = '<a class="pageref" href="%s">%s</a>' % (link, labelfmt %(figids[target]))
                ref_info.replace_self(raw(html, html, format='html'))
            else:
                ref_info.replace_self(Text(labelfmt % (figids[target])))


def clean_env(app):
    app.builder.env.i=1
    app.builder.env.figid_docname_map = {}

def setup(app):
    app.add_config_value('number_figures', True, True)
    app.add_config_value('figure_caption_prefix', "Figure", True)

    app.add_node(page_ref,
                 text=(skip_page_ref, None),
                 html=(skip_page_ref, None),
                 latex=(latex_visit_page_ref, None))

    app.add_role('page', XRefRole(nodeclass=page_ref))

    app.add_node(num_ref,
                 latex=(latex_visit_num_ref, None))

    app.add_role('num', XRefRole(nodeclass=num_ref))

    app.connect("builder-inited", clean_env)
    app.connect('doctree-read', doctree_read)
    app.connect('doctree-resolved', doctree_resolved)
