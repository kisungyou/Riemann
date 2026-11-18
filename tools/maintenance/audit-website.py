#!/usr/bin/env python3
"""Read-only audit of rendered Riemann website; emits JSON and exits nonzero on defects."""
import argparse
from collections import Counter
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
import xml.etree.ElementTree as ET

ap = argparse.ArgumentParser()
ap.add_argument('root', nargs='?', default=str(Path(__file__).resolve().parents[2] / 'docs'))
ap.add_argument('--version', help='Expected navigation version; defaults to the package DESCRIPTION.')
a = ap.parse_args()
root = Path(a.root).resolve()
if a.version:
    expected_version = a.version
else:
    description = (root.parent / 'DESCRIPTION').read_text()
    expected_version = re.search(r'^Version:\s*(\S+)', description, re.M).group(1)
errors, notes = [], []
mathjax_pages=[]; mathjax_urls=set()
class Page(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]; self.ids=set(); self.text=[]; self.nav_versions=[]; self.in_version=False
    def handle_starttag(self, tag, attrs):
        d=dict(attrs)
        if 'id' in d: self.ids.add(d['id'])
        if tag == 'small' and 'nav-text' in d.get('class', ''): self.in_version=True
        for key in ('href', 'src'):
            if key in d: self.links.append((d[key], tag, d))
    def handle_endtag(self, tag):
        if tag == 'small': self.in_version=False
    def handle_data(self, data):
        self.text.append(data)
        if self.in_version: self.nav_versions.append(data.strip())

pages={}
for f in sorted(root.rglob('*.html')):
    p=Page(); p.feed(f.read_text(errors='replace')); pages[f]=p

def resolve(source, url):
    u=urlsplit(url)
    if u.scheme or u.netloc:
        if u.scheme not in ('http', 'https') or u.netloc not in ('www.kisungyou.com','kisungyou.com'): return None
        if not u.path.startswith('/Riemann/'): return None
        target=root/unquote(u.path[len('/Riemann/'):])
    elif u.path.startswith('/'):
        if not u.path.startswith('/Riemann/'): return None
        target=root/unquote(u.path[len('/Riemann/'):])
    else:
        target=source.parent/unquote(u.path) if u.path else source
    target=target.resolve()
    if target.is_dir(): target=target/'index.html'
    return target, unquote(u.fragment)

def check(source, url, context):
    result=resolve(source,url)
    if result is None: return
    target,fragment=result
    if not target.exists(): errors.append({'type':'missing_target','source':str(source.relative_to(root)),'url':url,'context':context})
    elif fragment and target in pages and fragment not in pages[target].ids:
        errors.append({'type':'missing_anchor','source':str(source.relative_to(root)),'url':url,'context':context})

for f,p in pages.items():
    for url,tag,attrs in p.links:
        check(f,url,tag)
        mime=attrs.get('type','')
        if any(x in mime for x in ('“','”')):
            notes.append({'type':'curly_quoted_mime','source':str(f.relative_to(root)),'value':mime})
    if p.nav_versions and p.nav_versions != [expected_version]:
        errors.append({'type':'current_version','source':str(f.relative_to(root)),'value':p.nav_versions})
    content=' '.join(p.text)
    renderer=[url for url,tag,attrs in p.links if tag == 'script' and 'mathjax' in url.lower()]
    if renderer:
        mathjax_pages.append(str(f.relative_to(root))); mathjax_urls.update(renderer)
    if re.search(r'\\(?:\(|\[)',content) and not renderer:
        errors.append({'type':'tex_without_mathjax','source':str(f.relative_to(root))})
    if re.search(r'(?:/Users/|/private/|file://|/tmp/)',content):
        errors.append({'type':'private_path','source':str(f.relative_to(root))})
    if re.search(r'(?:stationary result.{0,20}local solution|stationarity.{0,20}certifies.{0,20}local minimum)',content,re.I):
        errors.append({'type':'stale_stationarity','source':str(f.relative_to(root))})

css_files=list(root.rglob('*.css'))
for css in css_files:
    for match in re.finditer(r'url\(\s*[\"\']?([^\"\'\)]+)',css.read_text(errors='replace')):
        check(css,match.group(1).strip(),'css_url')

sitemap=root/'sitemap.xml'
if not sitemap.exists(): errors.append({'type':'missing_sitemap'})
else:
    for el in ET.parse(sitemap).getroot().iter():
        if el.tag.endswith('loc') and el.text: check(sitemap,el.text,'sitemap')

search=root/'search.json'
if not search.exists(): errors.append({'type':'missing_search'})
else:
    entries=json.loads(search.read_text())
    for entry in entries:
        if not isinstance(entry,dict):continue
        url=entry.get('path')
        if isinstance(url,str):
            fragment=entry.get('id')
            check(search,url+('#'+fragment if isinstance(fragment,str) and fragment else ''),'search')
        text=entry.get('text','')
        if re.search(r'stationary result.{0,20}local solution',text,re.I):
            errors.append({'type':'stale_search_stationarity','url':url})
    mean_entries=[x for x in entries if isinstance(x,dict) and 'riem.mean.html' in str(x.get('path'))]
    if not any('saddle' in str(x).lower() and 'maximum' in str(x).lower() for x in mean_entries):
        errors.append({'type':'missing_corrected_stationarity_in_search'})

reference=root/'reference'/'riem.mean.html'
if reference not in pages or not all(w in ' '.join(pages[reference].text).lower() for w in ('saddle','maximum')):
    errors.append({'type':'missing_corrected_stationarity_reference'})

manifest=root/'site.webmanifest'
if manifest.exists():
    data=json.loads(manifest.read_text())
    for icon in data.get('icons',[]):
        if isinstance(icon,dict) and icon.get('src'):check(manifest,icon['src'],'manifest')

man=root.parent/'man'
missing_topics=[f.stem for f in man.glob('*.Rd') if not (root/'reference'/(f.stem+'.html')).exists()]
if missing_topics: errors.append({'type':'missing_reference_topics','topics':missing_topics})
print(json.dumps({'root':str(root),'html_pages':len(pages),'css_files':len(css_files),'rd_topics':len(list(man.glob('*.Rd'))),
 'mathjax_page_count':len(mathjax_pages),'mathjax_urls':sorted(mathjax_urls),
 'errors':errors,'error_count':len(errors),'informational_counts':dict(Counter(x['type'] for x in notes)),
 'informational_examples':notes[:3]},indent=2))
raise SystemExit(1 if errors else 0)
