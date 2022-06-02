#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#
#  Utilities for converting morph analyses from 
#  the UniversalDependencies (UD)  format  to a 
#  reduced Vabamorf's format
#
import os, os.path, re
from collections import defaultdict
from collections import OrderedDict

from estnltk.text import Text
from estnltk.layer.layer import Layer
from estnltk.layer.annotation import Annotation


# In[1]:



# =================================================
# =================================================
#    Convert UD annotations to reduced 
#    Vabamorf's annotations
# =================================================
# =================================================

def convert_ud_layer_to_reduced_morph_layer( text_obj, ud_layer, output_layer, add_layer=True ):
    '''Creates a reduced version of the UD layer which consists only of morph_analysis. 
       The reduced morph layer contains attributes 'lemma', 'pos' and 'form',
       and it uses Vabamorf's morphological categories.
    '''
    assert isinstance(text_obj, Text)
    #SEE OSA LOOB UUE TÜHJA KIHI, KUHU MÄRGENDUS LÄHEB
    assert ud_layer in text_obj.layers,            '(!) Layer {!r} missing from: {!r}'.format(ud_layer, text_obj.layers)
    redux_layer = Layer(name=output_layer,                   attributes=('lemma', 'pos', 'form'),                   text_object=text_obj,                  ambiguous=True)
    #SEE OSA PANEB MÄRGENDUSE SINNA KIHTI
    for ud_word in text_obj[ ud_layer ]:
        for ann in ud_word.annotations:
            attribs = parse_ud_morph_redux_attribs( ann )
            redux_layer.add_annotation( ud_word.base_span, **attribs )
            pass
    #SEE LISAB KIHI SÕNALE
    if add_layer:
        text_obj.add_layer( redux_layer )
    return redux_layer

def _split_feats( morph_form_feats ):
    '''Creates a dictionary based on UD's "feats" attribute.'''
    if morph_form_feats is None or len(morph_form_feats) == 0:
        return {}
    feat_chunks = morph_form_feats.split('|')
    feat_chunks_split = [chunk.split('=') for chunk in feat_chunks]
    feats = {kv[0]:kv[1] for kv in feat_chunks_split if len(kv) == 2}
    return feats


def _clean_lemma( lemma ):
    '''Removes '=' symbols from lemma if they appear between letters.'''
    new_lemma = []
    for i in range(len(lemma)):
        last_c = lemma[i-1] if i-1>-1 else ''
        c = lemma[i]
        next_c = lemma[i+1] if i+1<len(lemma) else ''
        if c == '=':
            if not(last_c.isalpha() and next_c.isalpha()):
                new_lemma.append( c )
        else:
            new_lemma.append( c )
    return ''.join(new_lemma)

# Mapping cases from UD to Vabamorf
ud_to_vm_case_mapping = {
    'Nom':'n', 
    'Gen':'g',
    'Par':'p',
    'Ill':'ill',
    'Ine':'in',
    'Ela':'el',
    'All':'all',
    'Ade':'ad',
    'Abl':'abl',
    'Tra':'tr',
    'Ter':'ter',
    'Ess':'es',
    'Abe':'ab',
    'Com':'kom',
    # aditiiv
    'Add':'adt'
}


def parse_ud_morph_redux_attribs( ud_annotation):
    '''Creates and returns a reduced version UD's morphosyntactic annotation.
       The reduced version contains only attributes 'lemma', 'pos' and 'form', 
       and uses Vabamorf's morphological categories.
       
       See also:
        https://github.com/EstSyntax/EstUD/blob/master/cgmorf2conllu/cgmorf2conllu.py
    '''
    attribs = { 'lemma':'', 'pos':'', 'form':'' }
    # UD attributes:
    ud_lemma = ud_annotation['lemma']
    #XPOSTAGID ON VABAMORFI TEHTUD ANNOTATSIOONID? UPOSTAGID ON UD ANNOTATSIOONID?
    #KUST XPOS JA UPOS TULI? SEDA POLNUD SEAL ud_syntax KIHIL LOETLETUD
    #MIKS MEIL ON VAJA TEADA, MIS VANA VABAMORFI TAG ON? 
    #KAS SEDA EI PEAKS LIHTSALT ISE VABAMORFIKS TEGEMA, MITTE SEDA USALDAMA?
    if 'xpos' in ud_annotation:
        ud_xpos  = ud_annotation['xpos'] # the old postag
    elif 'xpostag' in ud_annotation:
        ud_xpos  = ud_annotation['xpostag'] # the old postag
    if 'upos' in ud_annotation:
        ud_upos  = ud_annotation['upos'] # the new (UD) postag
    elif 'upostag' in ud_annotation:
        ud_upos  = ud_annotation['upostag'] # the new (UD) postag
    ud_feats = _split_feats(ud_annotation['feats']) if isinstance(ud_annotation['feats'], str) else ud_annotation['feats']
    assert isinstance(ud_feats, (dict, OrderedDict))
    # ==============================================================
    #   1) Parse lemma
    # ==============================================================
    attribs['lemma'] = ud_lemma
    if attribs['lemma'] is None:
        attribs['lemma'] = ''
    # Remove '=' symbols between letters:
    attribs['lemma'] = _clean_lemma( attribs['lemma'] )
    # ==============================================================
    #   2) Parse postag
    # ==============================================================
    attribs['pos']   = ud_xpos
    #siin muudame kõiki sõnaliike Vabamorfi vormile
    #siit on puudu VM genitiivatribuut(G) 
    ud_degree = ud_feats.get('Degree', None)
    ud_abbr = ud_feats.get('Abbr', None)
    
    if ud_degree == 'Pos':
        attribs['pos'] = 'A'
    if ud_degree == 'Cmp':
        attribs['pos'] = 'C'
    if ud_degree == 'Sup':
        attribs['pos'] = 'U'
    if ud_upos == 'ADV':
        attribs['pos'] = 'D'
    if ud_upos == 'PROPN':
        attribs['pos'] = 'H'
    if ud_upos == 'INTJ':
        attribs['pos'] = 'I'
    if ud_upos == 'CCONJ' or ud_upos == 'SCONJ':
        attribs['pos'] = 'J'
    if ud_upos == 'ADP':
        attribs['pos'] = 'K'
    if ud_upos == 'PRON' or ud_upos == 'DET':
        attribs['pos'] = 'P'
    if ud_upos == 'NOUN':
        attribs['pos'] = 'S'
    if ud_upos == 'VERB' or ud_upos == 'AUX':
        attribs['pos'] = 'V'
    if ud_upos == 'PUNCT' or ud_upos == 'SYM':
        attribs['pos'] = 'Z'
    if ud_abbr == 'Yes':
        attribs['pos'] = 'Y'
    if ud_upos == 'X':
        attribs['pos'] = 'X'
    if ud_xpos == 'G':
        attribs['pos'] = 'G'
    
    if ud_upos == 'NUM':
        ud_numtype = ud_feats.get('NumType', None)
        if ud_numtype == 'Card':
            attribs['pos'] = 'N'
        if ud_numtype == 'Ord':
            attribs['pos'] = 'O'
    elif ud_upos == 'ADJ':
        # Reduce adjectives to numbers iff required
        ud_numtype = ud_feats.get('NumType', None)
        if ud_numtype == 'Card':
            attribs['pos'] = 'N'
        if ud_numtype == 'Ord':
            attribs['pos'] = 'O'
    # Interjection:  B == I (actually: I is subtype of B) (specific to EWTB corpus)
    if ud_xpos == 'B':
        attribs['pos'] = 'I'
    # Emoticons:  E == Z (specific to EWTB corpus)
    if ud_xpos == 'E':
        attribs['pos'] = 'Z'
    
    # ==============================================================
    #   3) Parse form
    # ==============================================================
    #  Nominal: ud has both case and number
    if 'Number' in ud_feats and 'Case' in ud_feats:
        ud_number = ud_feats['Number']
        vm_number = 'pl' if ud_number == 'Plur' else ud_number
        vm_number = 'sg' if vm_number == 'Sing' else vm_number
        ud_case = ud_feats['Case']
        assert ud_case in ud_to_vm_case_mapping,                '(!) Unexpected case {!r} in: {!r}'.format(ud_case, ud_annotation)
        vm_case = ud_to_vm_case_mapping[ud_case]
        attribs['form'] = (vm_number+' '+vm_case) if vm_case != 'adt' else vm_case
        #MIKS SIIN SEE PROBLEEME TEKITAB? 
        #
        # Special case -- long illative -- leads to an ambiguity:
        #  ud_case == 'Ill' --> 'sg ill' or 'adt'
        #  ud_case == 'Add' --> 'sg ill' or 'adt'
        # TODO: do we need to generate several variants here?
        # 
    # ... All the dance with the verbs ...
    if ud_upos == 'VERB' or ud_upos == 'AUX':
        # Get UD's category values
        ud_verb_form   = ud_feats.get('VerbForm', None)     # Fin, Inf, Part, Sup, Conv
        ud_voice       = ud_feats.get('Voice',    None)     # Act, Pass
        ud_mood        = ud_feats.get('Mood',     None)     # Ind, Imp, Cnd, Qou
        ud_case        = ud_feats.get('Case',     None)     # Ill, Ine, Ela, Tra, Abe
        ud_number      = ud_feats.get('Number',   None)     # Plur, Sing
        ud_person      = ud_feats.get('Person',   None)     # 1, 2, 3
        ud_tense       = ud_feats.get('Tense',    None)     # Past, Pres
        ud_polarity    = ud_feats.get('Polarity', None)     # Neg
        ud_connegative = ud_feats.get('Connegative', None)  # Yes
        assert not (ud_upos == 'VERB' and ud_case != None and ud_number != None),                '(!) There should be no such verb: {!r}!'.format( ud_annotation )
        #
        #  For an overview of Vabamorf's verb categories, 
        #  see: http://www.filosoft.ee/html_morf_et/morfoutinfo.html#4
        #
        # V1) Infinite forms
        # pure infinite
        if ud_verb_form == 'Inf':
            attribs['form'] = 'da'
        # supine personal
        if ud_verb_form == 'Sup' and ud_voice == 'Act':
            if ud_case == 'Ill':
                attribs['form'] = 'ma'
            if ud_case == 'Ine':
                attribs['form'] = 'mas'
            if ud_case == 'Ela':
                attribs['form'] = 'mast'
            if ud_case == 'Tra':
                attribs['form'] = 'maks'
            if ud_case == 'Abe':
                attribs['form'] = 'mata'
        # supine impersonal
        if ud_verb_form == 'Sup' and ud_voice == 'Pass':
            attribs['form'] = 'tama'
        # nud/tud
        if ud_verb_form == 'Part' and ud_tense == 'Past':
            if ud_voice == 'Act':
                attribs['form'] = 'nud'
            if ud_voice == 'Pass' and ud_polarity == 'Neg':
                attribs['form'] = 'neg tud'
            if ud_voice == 'Pass':
                attribs['form'] = 'tud'
        # v/tav
        if ud_verb_form == 'Part' and ud_tense == 'Pres':
            if ud_voice == 'Act':
                attribs['form'] = 'v'
            if ud_voice == 'Pass':
                attribs['form'] = 'tav'
        
        # ger
        if ud_verb_form == 'Conv':
            attribs['form'] = 'des'
        # V2) Negatives:
        if ud_polarity == 'Neg' or ud_connegative == 'Yes':
            # neg auxiliary
            if ud_upos == 'AUX' and ud_lemma in ['ära', 'ei']:
                attribs['form'] = 'neg'
            # neg personal 
            if ud_voice == 'Act':
                # Ind, Imp, Cnd, Qou
                if ud_mood == 'Ind' and ud_tense == 'Pres':
                    # (!) Ambiguity:  vm_form in ['o', 'neg o']
                    attribs['form'] = 'neg o'
                if ud_mood == 'Ind' and ud_tense == 'Past':
                    # (!) Ambiguity:  vm_form in ['nud', 'neg nud']
                    attribs['form'] = 'neg nud'
                #LIHTSALT O MÄRGEND ON JAATAVA KÕNE JAOKS, MIKS SEE SIIN ON?
                #if ud_mood == 'Imp' and ud_tense == 'Pres' and ud_person == '2' and ud_number == 'Sing':
                    #attribs['form'] = 'o'
                if ud_mood == 'Imp' and ud_tense == 'Pres' and ud_person == '2' and ud_number == 'Plur':
                    attribs['form'] = 'neg ge' 
                if ud_mood == 'Imp' and ud_tense == 'Pres' and ud_person == '3' and ud_number == 'Plur':
                    attribs['form'] = 'neg gu'
                #NEG GEM (ÄRGEM) ON TÄPSELT SAMA MIS NEG ME (ÄRME)
                #if ud_mood == 'Imp' and ud_tense == 'Pres' and ud_person == '1' and ud_number == 'Plur':
                    #attribs['form'] = 'neg me'
                if ud_mood == 'Imp' and ud_tense == 'Pres' and ud_person == '1' and ud_number == 'Plur':
                    attribs['form'] = 'neg gem'
                if ud_mood == 'Cnd' and ud_tense == 'Pres':
                    # (!) Ambiguity:  vm_form in ['ks', 'neg ks']
                    attribs['form'] = 'neg ks'
                if ud_mood == 'Cnd' and ud_tense == 'Past':
                    attribs['form'] = 'neg nuks'
                if ud_mood == 'Qot' and ud_tense == 'Pres':
                    attribs['form'] = 'neg vat'
            # neg impersonal 
            if ud_voice == 'Pass':
                if ud_mood == 'Imp' and ud_tense == 'Pres':
                    attribs['form'] = 'neg gu'
                if ud_mood == 'Ind' and ud_tense == 'Pres':
                    attribs['form'] = 'ta'
                
                
        ud_affirmative = (not ud_polarity == 'Neg') and (not ud_connegative == 'Yes')
        # V3) Indicative, affirmative
        if ud_affirmative and ud_mood == 'Ind':
            # Present tense
            if ud_number == 'Sing'   and ud_tense == 'Pres' and ud_person == '1':
                attribs['form'] = 'n'
            if ud_number == 'Plur'   and ud_tense == 'Pres' and ud_person == '1':
                attribs['form'] = 'me'
            if ud_number == 'Sing'   and ud_tense == 'Pres' and ud_person == '2':
                attribs['form'] = 'd'
            if ud_number == 'Plur'   and ud_tense == 'Pres' and ud_person == '2':
                attribs['form'] = 'te'
            if ud_number == 'Sing'   and ud_tense == 'Pres' and ud_person == '3':
                attribs['form'] = 'b'
            if ud_number == 'Plur'   and ud_tense == 'Pres' and ud_person == '3':
                attribs['form'] = 'vad'
            # Passive voice
            if ud_voice == 'Pass' and ud_tense == 'Pres' and ud_person == None:
                attribs['form'] = 'takse'
            # Past tense
            if ud_number == 'Sing'  and ud_tense == 'Past' and ud_person == '1':
                attribs['form'] = 'sin'
            if ud_number == 'Plur'  and ud_tense == 'Past' and ud_person == '1':
                attribs['form'] = 'sime'
            if ud_number == 'Sing'  and ud_tense == 'Past' and ud_person == '2':
                attribs['form'] = 'sid'
            if ud_number == 'Plur'  and ud_tense == 'Past' and ud_person == '2':
                attribs['form'] = 'site'
            if ud_number == 'Sing'  and ud_tense == 'Past' and ud_person == '3':
                attribs['form'] = 's'
            if ud_number == 'Plur'  and ud_tense == 'Past' and ud_person == '3':
                attribs['form'] = 'sid'
            # Passive voice
            if ud_voice == 'Pass' and ud_tense == 'Past':
                attribs['form'] = 'ti'
        # V4) Imperative, affirmative
        if ud_affirmative and ud_mood == 'Imp':
            #GU ON JU AINULT KOLMANDA ISIKU JAOKS? MIKS ON ERALDI PERSON==NONE, aga singular?
            #Active voice
            if ud_number == 'Sing'  and ud_tense == 'Pres' and ud_person == None and ud_voice == 'Act':
                attribs['form'] = 'gu'
            if ud_number == 'Sing'  and ud_tense == 'Pres' and ud_person == '2' and ud_voice == 'Act':
                attribs['form'] = 'o'
            if ud_number == 'Sing'  and ud_tense == 'Pres' and ud_person == '3' and ud_voice == 'Act':
                attribs['form'] = 'gu'
            if ud_number == 'Plur'  and ud_tense == 'Pres' and ud_person == '1' and ud_voice == 'Act':
                attribs['form'] = 'gem'
            if ud_number == 'Plur'  and ud_tense == 'Pres' and ud_person == '2' and ud_voice == 'Act':
                attribs['form'] = 'ge'
            if ud_number == 'Plur'  and ud_tense == 'Pres' and ud_person == '3' and ud_voice == 'Act':
                attribs['form'] = 'gu'
            #Passive voice
            if ud_voice == 'Pass' and ud_tense == 'Pres':
                attribs['form'] = 'tagu'
        # V5) Quotative, affirmative
        if ud_affirmative and ud_mood == 'Qot':
            if ud_tense == 'Pres' and ud_voice == 'Act':
                attribs['form'] = 'vat'
            if ud_tense == 'Pres' and ud_voice == 'Pass':
                attribs['form'] = 'tavat'
            if ud_tense == 'Past' and ud_voice == 'Act' :
                attribs['form'] = 'nuvat'
            if ud_tense == 'Past' and ud_voice == 'Pass':
                attribs['form'] = 'tuvat'
        # V6) Conditional, affirmative
        if ud_affirmative and ud_mood == 'Cnd':
            # Present tense
            if ud_tense == 'Pres' and ud_voice == 'Act' and ud_number == 'Sing' and ud_person == '1':
                # (!) Ambiguity:  vm_form in ['ksin', 'ks']
                attribs['form'] = 'ksin'
            if ud_tense == 'Pres' and ud_voice == 'Act' and ud_number == 'Sing' and ud_person == '2':
                # (!) Ambiguity:  vm_form in ['ksid', 'ks']
                attribs['form'] = 'ksid'
            if ud_tense == 'Pres' and ud_voice == 'Act' and ud_number == 'Sing' and ud_person == '3':
                attribs['form'] = 'ks'
            if ud_tense == 'Pres' and ud_voice == 'Act' and ud_number == 'Plur' and ud_person == '1':
                # (!) Ambiguity:  vm_form in ['ksime', 'ks']
                attribs['form'] = 'ksime'
            if ud_tense == 'Pres' and ud_voice == 'Act' and ud_number == 'Plur' and ud_person == '2':
                # (!) Ambiguity:  vm_form in ['ksite', 'ks']
                attribs['form'] = 'ksite'
            if ud_tense == 'Pres' and ud_voice == 'Act' and ud_number == 'Plur' and ud_person == '3':
                # (!) Ambiguity:  vm_form in ['ksid', 'ks']
                attribs['form'] = 'ksid'
            #EI OLE SELLIST KS-I, MILLE PERSON==NONE JU?
            if ud_voice == 'Act'  and ud_tense == 'Pres' and ud_person == None:
                attribs['form'] = 'ks'

            # Past tense
            #SAMA KÜSIMUS, MIKS ME PIKKA VORMI EELISTAME?
            if ud_tense == 'Past' and ud_voice == 'Act' and ud_number == 'Sing' and ud_person == '1':
                # (!) Ambiguity:  vm_form in ['nuksin', 'nuks']
                attribs['form'] = 'nuksin'
            if ud_tense == 'Past' and ud_voice == 'Act' and ud_number == 'Sing' and ud_person == '2':
                # (!) Ambiguity:  vm_form in ['nuksid', 'nuks']
                attribs['form'] = 'nuksid'
            if ud_tense == 'Past' and ud_voice == 'Act' and ud_number == 'Sing' and ud_person == '3':
                attribs['form'] = 'nuks'
            if ud_tense == 'Past' and ud_voice == 'Act' and ud_number == 'Plur' and ud_person == '1':
                # (!) Ambiguity:  vm_form in ['nuksime', 'nuks']
                attribs['form'] = 'nuksime'
            if ud_tense == 'Past' and ud_voice == 'Act' and ud_number == 'Plur' and ud_person == '2':
                # (!) Ambiguity:  vm_form in ['nuksite', 'nuks']
                attribs['form'] = 'nuksite'
            if ud_tense == 'Past' and ud_voice == 'Act' and ud_number == 'Plur' and ud_person == '3':
                # (!) Ambiguity:  vm_form in ['nuksid', 'nuks']
                attribs['form'] = 'nuksid'
            
            #Passive voice
            if ud_tense == 'Pres' and ud_voice == 'Pass':
                attribs['form'] = 'taks'
            if ud_tense == 'Past' and ud_voice == 'Pass':
                attribs['form'] = 'tuks'
        
    return attribs

