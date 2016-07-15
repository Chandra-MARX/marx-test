.. _sect-{{ module.__name__ }}:

{{ module.title|default(module.__name__) }}
=============================================================================
{{ module.__doc__ }}

{% for testinstance in testinstances %}

.. _sect-test-{{ testinstance.name }}:

{{ testinstance.title|default(testinstance.name) }}
------------------------------------------------------------------------------

{% if testinstance.obsid %}
using data from: `ObsID {{ testinstance.obsid }} <http://cda.harvard.edu/chaser/startViewer.do?menuItem=details&obsid={{testinstance.obsid }}>`_
{% endif %}
link to code for this test: :ref:`test-code-{{testinstance.name}}`

{{ testinstance.__doc__ }}

{% for figname, figure in testinstance.figures.iteritems() %}

.. figure:: {{ figpath }}/{{testinstance.name}}_{{ figname }}.*
   :align: center
   :alt: {{ figure.alternative }}

   {{ figure.caption }}

{% endfor %}

{{ testinstance.summary }}

{% endfor %}
